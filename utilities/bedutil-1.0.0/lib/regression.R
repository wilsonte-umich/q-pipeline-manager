
#######################################################################
# 'regression.R' takes a BED file from STDIN to which a second score
# has been appended as the 7th field.  A correlation plot is made and 
# a linear regression performed.  The original data are repeated to
# STDOUT and the regression data are printed to STDERR.
#----------------------------------------------------------------------
# Usage:
#     Rscript regression.R <jpgFile> <xLabel> <yLabel> <maxInputScore> <maxCalcScore>
# where the required parameters are:
#     jpgFile        name of jpg graph to be created
#     xLabel         label for the X-axis, i.e. the name of the input scores
#     yLabel         label for the Y-axis, i.e. the name of the appended scores
#     maxInputScore  maximum allowed value for the input scores
#                    any higher scores are set to maxInputScore prior to regression
#                    a value of zero enforces no limit
#     maxCalcScore   maximum allowed value for the appended (calculated) scores
#                    any higher scores are set to maxCalcScore prior to regression
#                    a value of zero enforces no limit
#######################################################################

# extract passed parameters
args <- commandArgs(TRUE)
jpgFile <- args[1]
xLabel <- args[2]
yLabel <- args[3]
maxInputScore <- as.numeric(args[4])
maxCalcScore <- as.numeric(args[5])

# collect and parse the input data
input <- read.table("stdin",header=FALSE,sep="\t")
inputScores <- input[,'V5']
calcScores <- input[,'V7']
if(maxInputScore != 0){
    inputScores <- sapply(inputScores, min, maxInputScore)
}
if(maxCalcScore != 0){
    calcScores <- sapply(calcScores, min, maxCalcScore) 
}
maxAdjInputScore <- max(inputScores)
maxAdjCalcScore <- max(calcScores)
solidPoints <- input[input$V5<=maxAdjInputScore & input$V7<=maxAdjCalcScore,c('V5','V7')]
openPoints  <- input[input$V5>maxAdjInputScore | input$V7>maxAdjCalcScore,c('V5','V7')]
if(maxInputScore != 0){
    openPoints$V5 <- sapply(openPoints$V5, min, maxInputScore)
}
if(maxCalcScore != 0){
    openPoints$V7 <- sapply(openPoints$V7, min, maxCalcScore)
}

# perform the linear regression
reg <- lm(calcScores~inputScores)

# create the visualization plot
width <- 2.5
height <- 2.5
units <- 'in' #w and h in inches
pointsize <- 7
res <- 600 #dpi
pch <- 20
cex <- 0.3
xlim <- c(0,maxAdjInputScore)
ylim <- c(0,maxAdjCalcScore)
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(solidPoints[,'V5'],solidPoints[,'V7'],xlab=xLabel,ylab=yLabel,pch=20,xlim=xlim,ylim=ylim)
points(openPoints[,'V5'],openPoints[,'V7'],pch="+")
abline(reg=reg)

# print the results
write.table(input,file="",quote=FALSE,sep="\t",na="",row.names=FALSE,col.names=FALSE)
sink(stderr())
summary(reg)

