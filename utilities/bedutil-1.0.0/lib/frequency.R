
#######################################################################
# 'simulate.R' takes tab-delimited score-frequency pairs from STDIN,
# fits a Gaussian to them, and creates a visualization plot.  The input
# data must be provided in the following format, including the header:
#     aggregateScore<TAB>frequency<TAB>fractionAtOrBelow<TAB>fractionAtOrAbove
#     <VALUE><TAB>-2<TAB>-2<TAB>-2                 (the offset control value, if any)
#     <VALUE><TAB>-1<TAB><VALUE><TAB><VALUE>       (the actual value)
#     <VALUE><TAB><VALUE><TAB><VALUE><TAB><VALUE>  (frequency-ordered simulation distribution)
# All input lines are repeated to STDOUT.  The following simulation values
# are printed to STDERR in the indicated format:
#     <VALUE><TAB>amplitude
#     <VALUE><TAB>mean
#     <VALUE><TAB>stddev (the first three values relate the simulation curve-fit)
#     <VALUE><TAB>p-value (2-sided p-value of the actual aggregateScore, from the curve-fit)
#     <VALUE><TAB>freq-tail (1-sided frequency of iterations at or beyond the actual aggregateScore)
#----------------------------------------------------------------------
# Usage:
#     Rscript simulate.R <jpgFile> <xLabel>
# where:
#     jpgFile    name of jpg graph to be created
#     xLabel     label for the X-axis, i.e. the name of the score
#######################################################################

# extract passed parameters
args <- commandArgs(TRUE)
jpgFile <- args[1]
xlab <- args[2]

# collect and parse the input data
input <- read.table("stdin",header=TRUE,sep="\t")
write.table(input,file="",quote=FALSE,sep="\t",na="",row.names=FALSE,col.names=TRUE)  # repeat to STDOUT
offsetControl <- input[input$frequency==-2,'aggregateScore']
actual <- input[input$frequency==-1,c('aggregateScore','fractionAtOrBelow','fractionAtOrAbove')]
data <- input[input$frequency>=0,c('aggregateScore','frequency')]
maxFreq <- max(data$frequency,na.rm=TRUE)
a0 <- maxFreq
m0 <- sum(data$frequency * data$aggregateScore) 
data$freqSum <- cumsum(data$frequency)
s0 <- min(data[data$freqSum>=0.841,'aggregateScore']) - m0
noOffset <- -9999.9999
offsetControl <- ifelse(length(offsetControl)==0, noOffset, offsetControl[1])

# try to fit a Gaussian curve to the simulation frequencies
guassian <- frequency ~ a*(exp(-((aggregateScore-m)^2)/(2*(s^2))))
curveFit <- try(nls(guassian,data=data,start=list(a=a0,m=m0,s=s0)), silent=TRUE)
curveFitOK <- !inherits(curveFit, "try-error")
if(curveFitOK){
    coef <- coef(curveFit)
    a <- coef[['a']]
    m <- coef[['m']]
    s <- coef[['s']]
} else {
    a <- a0
    m <- m0
    s <- s0
}

# determine the 2-sided p-value of the actual value, as well as
# the 1-sided frequency of values at or beyond the actual value
actualAgg <- actual[1,'aggregateScore']
if(actualAgg<=m){
    pval <- ifelse(curveFitOK, pnorm(actualAgg,mean=m,sd=s,lower.tail=TRUE) * 2, "na")
    tail <- actual[1,'fractionAtOrBelow']
} else {
    pval <- ifelse(curveFitOK, pnorm(actualAgg,mean=m,sd=s,lower.tail=FALSE) * 2, "na")
    tail <- actual[1,'fractionAtOrAbove']
}

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
main <- paste("p=", ifelse(pval=="na", pval, sprintf("%.1E",pval)), ", f=", ifelse(tail==0,0,sprintf("%.1E",tail)), sep="")

# set axis limits
f05 <- data[data$frequency>=maxFreq*0.05,]
offsetWorking <- ifelse(offsetControl==noOffset, actualAgg, offsetControl)
xmin <- min(m-3.5*s, actualAgg, offsetWorking, min(f05$aggregateScore))/1.1
xmax <- max(m+3.5*s, actualAgg, offsetWorking, max(f05$aggregateScore))*1.1
xlim <- c(xmin,xmax)
ymax <- maxFreq * 1.2
ylim <- c(0,ymax)

# plot the data and fit
for (ext in c("",".pdf")){ # create both jpg and a pdf image
    imageFile <- paste(jpgFile, ext, sep="")
    if(ext==""){
        jpeg(file=imageFile,width=width,height=height,unit=units,pointsize=pointsize,res=res)
    }else{
        pdf(file=imageFile,width=width,height=height,pointsize=pointsize)
    }
    plot(data$aggregateScore,data$frequency,type='h',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col='blue')
    title(main,col.main="darkgreen",font.main=1)
    if(curveFitOK){
        x <- seq(from=xmin,to=xmax,length.out=1000)
        y <- a*(exp(-((x-m)^2)/(2*(s^2))))
        lines(x,y,col='red')
    }
    lines(c(actualAgg,actualAgg),ylim,col="darkgreen")
    if(offsetControl!=noOffset){ lines(c(offsetControl,offsetControl),ylim) }
    graphics.off()
}

# print the results
sink(stderr())
cat(a0,   "\t",'amplitude-guess',"\n",sep="")
cat(m0,   "\t",'mean-guess',     "\n",sep="")
cat(s0,   "\t",'stddev-guess',   "\n",sep="")
cat(a,   "\t",'amplitude',"\n",sep="")
cat(m,   "\t",'mean',     "\n",sep="")
cat(s,   "\t",'stddev',   "\n",sep="")
cat(pval,"\t",'p-value',  "\n",sep="")
cat(tail,"\t",'frequency-tail',  "\n",sep="")

