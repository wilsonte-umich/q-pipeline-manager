
#######################################################################
# 'relate' carries out the work of relate.pl, i.e. 'bedutil relate'
#######################################################################

# functions
asLogical <- function(boolean){ return( boolean=="1" | boolean=="TRUE" ) }
getEnv <- function(name, default){
    v<-Sys.getenv(name)
    return(ifelse(v=="",default,v))
}
applyLog <- function(x, isLog){
    return( ifelse(isLog, ifelse(x<=0,NA,log10(x)), x) )
}
purgeNA <- function(df){
    df <- df[!is.na(df[[cx]])&!is.na(df[[cy]]),]
    return(df)
}

# extract passed parameters
cat("setting parameters\n")
JPG_FILE <- Sys.getenv("JPG_FILE")
INDEPENDENT <- Sys.getenv("INDEPENDENT")
DEPENDENT <- Sys.getenv("DEPENDENT")
SUPPRESS_SIM <- asLogical(Sys.getenv("SUPPRESS_SIM"))
X_LOG <- asLogical(Sys.getenv("X_LOG"))
X_MIN <- Sys.getenv("X_MIN")
X_MAX <- Sys.getenv("X_MAX")
Y_LOG <- asLogical(Sys.getenv("Y_LOG"))
Y_MIN <- Sys.getenv("Y_MIN")
Y_MAX <- Sys.getenv("Y_MAX")
ROUND_DIGITS <- Sys.getenv("ROUND_DIGITS")

# set common plot properties
width <- 3.5  # could add getEnv to allow any value to be passed via environment
height <- 3.167
units <- 'in' #w and h in inches
pointsize <- 8
res <- 600 #dpi
quality <- 50
actualPointColor <- "grey"
actualColor <-      "red"
simulationColor <-  "darkgreen"
simulationColor2 <- "blue"
meanLty <- 1
medianLty <- 2
quantLty <- 3
actualPch <- as.numeric(getEnv("ACTUAL_PCH",1))
actualCex <- as.numeric(getEnv("ACTUAL_CEX",1))

# collect and parse the input data
cat("reading data\n")
data <- read.table("stdin",header=FALSE,sep="\t",stringsAsFactors=FALSE)
cx <- 5
cy <- ncol(data)
ci <- cy - 1
d <-  data[,c(ci,cx,cy)]
ci <- 1
cx <- 2
cy <- 3

# purge non-numeric data for some scores
cat("removing NULL scores\n")
d <- d[d[[cy]]!='NULL',]
d$ny <- as.numeric(d[[cy]])
cy <- 4

# calculate the correlation coefficients on logged but not-limited X and Y data
cat("calculating correlation coefficients\n")
actual <- d[d[[ci]]==-1,]
actual[[cx]] <- sapply(actual[[cx]],applyLog,X_LOG)
actual[[cy]] <- sapply(actual[[cy]],applyLog,Y_LOG)
actual <- purgeNA(actual)
pearson <-  round( cor(actual[[cx]], actual[[cy]], method="pearson")  ,3)
spearman <- round( cor(actual[[cx]], actual[[cy]], method="spearman") ,3)
main <- paste("r(p)=", pearson, ", r(s)=", spearman, sep="")
cat(paste(pearson,  "\tpearson correlation coefficient\n",  sep=""))
cat(paste(spearman, "\tspearman correlation coefficient\n", sep=""))

# adjust x values so that all points appear on graph
# outliers are aggregated into min and max bins
cat("applying x limits\n")
xmin <- ifelse(X_LOG,min(d[d[[cx]]>0,cx],na.rm=TRUE),0)
xmin <- ifelse(X_MIN=="",xmin,as.numeric(X_MIN))
xmax <- ifelse(X_MAX=="",max(d[[cx]],na.rm=TRUE),as.numeric(X_MAX))
d[[cx]] <- ifelse(d[[cx]]<xmin,xmin,d[[cx]])
d[[cx]] <- ifelse(d[[cx]]>xmax,xmax,d[[cx]])

# log the data if instructed
cat("applying logs (if any)\n")
d[[cx]] <- sapply(d[[cx]],applyLog,X_LOG)
xmin <- applyLog(xmin,X_LOG)
xmax <- applyLog(xmax,X_LOG)
d[[cy]] <- sapply(d[[cy]],applyLog,Y_LOG)
d <- purgeNA(d)

# create the x bins for aggregation
cat("binning x data\n")
if(ROUND_DIGITS==""){
    d$xb <- d[[cx]]
}else{
    d$xb <- round(d[[cx]],digits=as.numeric(ROUND_DIGITS))
}
cxb <- 5
xmin <- min(xmin, d[[cxb]], na.rm=TRUE) # in case bin rounding increased the x limits
xmax <- max(xmax, d[[cxb]], na.rm=TRUE)

# separate actual from simulated data
actual <- d[d[[ci]]==-1,]
if(!SUPPRESS_SIM){ simulation <- d[d[[ci]]>=0,] }

# aggregate the separated data
# aggregates calculated using logged but not-limited Y data, and logged and limited X data
cat("calculating aggregates\n")
aY <-        actual[[cy]]
aBins <-     list(bin=actual[[cxb]])
actMean <-   aggregate(aY, aBins, mean,   na.rm=TRUE)
actMedian <- aggregate(aY, aBins, median, na.rm=TRUE)
if(!SUPPRESS_SIM){
    sY <-          simulation[[cy]]
    sBins <-       list(bin=simulation[[cxb]])
    # aggregate all iteration values within the bin, i.e. a straight mean/median
    simMean <-     aggregate(sY, sBins, mean,   na.rm=TRUE)
    simMedian <-   aggregate(sY, sBins, median, na.rm=TRUE)
    simQ5 <-       aggregate(sY, sBins, quantile, probs=0.05, na.rm=TRUE)
    simQ95 <-      aggregate(sY, sBins, quantile, probs=0.95, na.rm=TRUE)
    # aggregate as the mean of mean, median of medians
    sItBins <-     list(iter=simulation[[ci]],bin=simulation[[cxb]])
    simItMean <-   aggregate(sY, sItBins, mean,   na.rm=TRUE)
    simItMedian <- aggregate(sY, sItBins, median, na.rm=TRUE)
    meanY <-       simItMean[,3]
    medianY <-     simItMedian[,3]
    meanBins <-    list(bin=simItMean$bin)
    medianBins <-  list(bin=simItMedian$bin)
    simMeanMean <-      aggregate(meanY,   meanBins,   mean,   na.rm=TRUE)
    simMedianMedian <-  aggregate(medianY, medianBins, median, na.rm=TRUE)
    simMedianQ5 <-      aggregate(medianY, medianBins, quantile, probs=0.05, na.rm=TRUE)
    simMedianQ95 <-     aggregate(medianY, medianBins, quantile, probs=0.95, na.rm=TRUE)
}

# adjust y values using provided limits
# outliers will plot at limits (but were aggregated using actual values)
cat("applying y limits\n")
ymin <- ifelse(Y_LOG,min(d[d[[cy]]>0,cy],na.rm=TRUE),0)
ymin <- ifelse(Y_MIN=="",ymin,as.numeric(Y_MIN))
ymax <- ifelse(Y_MAX=="",max(actual[[cy]],na.rm=TRUE),as.numeric(Y_MAX))
actual[[cy]] <- ifelse(actual[[cy]]<ymin,ymin,actual[[cy]])
actual[[cy]] <- ifelse(actual[[cy]]>ymax,ymax,actual[[cy]])

# set axes
cat("generating plots\n")
xlab <- INDEPENDENT
xlim <- c(xmin,xmax)
ylab <- DEPENDENT
ylim <- c(ymin,ymax)

# create both jpg and a pdf image
for (ext in c("",".pdf")){

# initialize plot image
imageFile <- paste(JPG_FILE, ext, sep="")
if(ext==""){
    jpeg(filename=imageFile,res=res,width=width,height=height,pointsize=pointsize,units=units)
}else{
    pdf(file=imageFile,width=width,height=height,pointsize=pointsize)
}
plot(xlim,ylim,type="n",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main)

# plot the actual points
cxy <- c(cx,cy)
points(actual[,cxy],pch=actualPch,cex=actualCex,col=actualPointColor)

# plot the aggregate and range lines
cxy <- c(1,2)
if(!SUPPRESS_SIM){
    lines(simQ5[,cxy],    col=simulationColor,lty=quantLty)
    lines(simQ95[,cxy],   col=simulationColor,lty=quantLty)
    lines(simMeanMean[,cxy],    col=simulationColor2,lty=meanLty)
    lines(simMedianMedian[,cxy],col=simulationColor2,lty=medianLty)
    lines(simMedianQ5[,cxy],    col=simulationColor2,lty=quantLty)
    lines(simMedianQ95[,cxy],   col=simulationColor2,lty=quantLty)
}
lines(actMean[,cxy],  col=actualColor,lty=meanLty)
lines(actMedian[,cxy],col=actualColor,lty=medianLty)

}

