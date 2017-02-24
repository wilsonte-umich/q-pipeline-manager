
#######################################################################
# 'distribution.R' does the work of command 'distribution'
#######################################################################

# extract passed parameters
IMAGE_PREFIX <- Sys.getenv("IMAGE_PREFIX")
X_LABEL <- Sys.getenv("X_LABEL")
SUPPRESS_NULL <- Sys.getenv("SUPPRESS_NULL")
N_BINS <- as.numeric(Sys.getenv("N_BINS"))
AGGREGATE_FUNCTION <- Sys.getenv("AGGREGATE_FUNCTION")

# collect and parse the input data (as formatted by crosstab)
input <- read.table("stdin",header=TRUE,sep="\t",stringsAsFactors=FALSE,na.strings=c("NA","NULL",""))

# handle null suppression
if(SUPPRESS_NULL==""){ input[is.na(input)] <- 0 }

# split actual and simulated data
as <- input[,c(7:ncol(input))]
a <- input[,7]
s <- input[,c(8:ncol(input))]

# aggregate the simulated data
t <- apply(s, 1, AGGREGATE_FUNCTION, na.rm=TRUE)

# calculate the paired differences
d  <- a - s
dt <- a - t

# percentile the actualdata using the simulation data
percentile <- function(v){ va<-v[[1]];if(is.na(va)){return(NA)};vs<-v[2:length(v)];vs<-vs[!is.na(vs)];vsl<-length(vs);if(vsl==0){return(NA)}; round((length(vs[vs<va]) + length(vs[vs==va])/2) / vsl * 100) }
p <-apply(as, 1, percentile)

# simplify simulation values to a vector
s <- c(as.matrix(s))
d <- c(as.matrix(d))

# initialize the visualization plot
ylab <- 'Frequency'
width <- 3.5
height <- 3.167
units <- 'in' #w and h in inches
pointsize <- 8
res <- 300 #dpi
pch <- 20
cex <- 0.3
quality <- 50

# set plot properties
col1 <- "red"
col2 <- "blue"
legendLwd <- 1.5
legendCex <- 0.8

# plotting function
makePlots <- function(d1,d2,is2,paired){

# set more plot properties
xlab <- ifelse( is2, X_LABEL, paste("delta(",X_LABEL,")",sep="") )

# bin the scores
xmin <- min(d1,na.rm=TRUE)
xmax <- max(d1,na.rm=TRUE)
if(is2){
    xmin <- min(xmin,d2,na.rm=TRUE)
    xmax <- max(xmax,d2,na.rm=TRUE)
}
binSize <- (xmax - xmin) / N_BINS
if(binSize==0){ binSize <- 1 }
b1 <- round(d1/binSize) * binSize
if(is2){ b2 <- round(d2/binSize) * binSize }
xmin <- min(b1,na.rm=TRUE)
xmax <- max(b1,na.rm=TRUE)
if(is2){
    xmin <- min(xmin,b2,na.rm=TRUE)
    xmax <- max(xmax,b2,na.rm=TRUE)
}

# calculate frequency histograms
cx <- 1
cy <- 2
h1 <- aggregate(b1,list(bin=b1),length)
h1[[cy]] <- h1[[cy]] / length(b1)
if(is2){
    h2 <- aggregate(b2,list(bin=b2),length)
    h2[[cy]] <- h2[[cy]] / length(b2)
}

# set up bar interleave
if(is2){
    quarterBinSize <- binSize / 4
    h1[[cx]] <- h1[[cx]] - quarterBinSize
    h2[[cx]] <- h2[[cx]] + quarterBinSize
    xmin <- xmin - quarterBinSize
    xmax <- xmax + quarterBinSize
}

# set axis limits
xlim <- c(xmin,xmax)
ymin <- 0
ymax <- max(h1[[cy]])
if(is2){ ymax <- max(ymax, h2[[cy]]) }
bpat <- 1.2*ymax
bpat2 <- c(1.067,1.133)*ymax
boxwex <- 0.11*ymax
boxwex2 <- 0.055*ymax
ymax <- ymax * 1.4
ylim <- c(ymin,ymax)
barsLwd <- 100 / ( (xmax-xmin) / binSize )

# plot the results
for (ext in c("jpg","pdf")){ # create both jpg and a pdf image

    # initalize the image
    imageFile <- paste(IMAGE_PREFIX, paired, ext, sep=".")
    if(ext=="jpg"){
        jpeg(file=imageFile,width=width,height=height,unit=units,pointsize=pointsize,res=res)
    }else{
        pdf(file=imageFile,width=width,height=height,pointsize=pointsize)
    }
    plot(xlim,ylim,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)

    # add a zero line to paired plots
    if(!is2){ lines(c(0,0),ylim) }

    # plot the histogram
    points(h1[[cx]],h1[[cy]],type='h',lwd=barsLwd,col=col1)
    if(is2){ points(h2[[cx]],h2[[cy]],type='h',lwd=barsLwd,col=col2) }

    # add the box plots
    if(is2){
        boxplot(d1,d2,border=c(col1,col2),horizontal=TRUE,add=TRUE,outline=FALSE,at=bpat2,boxwex=boxwex2)
    } else {
        boxplot(d1,border=col1,horizontal=TRUE,add=TRUE,outline=FALSE,at=bpat,boxwex=boxwex)
    }

    # add a legend
    if(is2){
        legend("topleft",col=c(col1,col2),lty=1,lwd=legendLwd,legend=c("actual","simulated"),bg="white",cex=legendCex)
    } else {
        legend("topleft",col=col1,lty=1,lwd=legendLwd,legend=c("actual - simulated"),bg="white",cex=legendCex)
    }

    graphics.off()
}

}

# create the plots
makePlots(a,s,TRUE,"unpaired_all")
makePlots(a,t,TRUE,"unpaired_agg")
makePlots(d,FALSE,FALSE,"paired_all")
makePlots(dt,FALSE,FALSE,"paired_agg")
makePlots(p,FALSE,FALSE,"percentile")
