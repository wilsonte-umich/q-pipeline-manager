
#######################################################################
# 'normtest.R' does the work of command 'normtest'
#######################################################################

# collect and parse the input data
source(paste(Sys.getenv("SCRIPT_DIR"), "/common.R", sep=""))
loadCrosstab(TRUE)  # sorted by size
checkIterations()

# perform normality test for all input features
maxnv <- 1000  # shapiro.test takes at most 5000 points
normtest <- function(row, log){
    feat <- getFeature(row, log)
    if(feat$nv<3){ return(0) }
    if(feat$nv>maxnv){ 
        feat$vs <- sort(sample(feat$vs, maxnv)) 
        feat$nv <- maxnv
    }     
    if (feat$vs[1] == feat$vs[feat$nv]){ return(0) } 
    sh <- shapiro.test(feat$vs)
    return(sh$p.value)
}
ips <-     sapply(1:nF, normtest, FALSE)
ips_log <- sapply(1:nF, normtest, TRUE)

# initialize the plot
main <- env['SIM_LABEL']
xlab <- "Feature Length"
xmin <- min(lens)
xmax <- max(lens)
xlim <- c(xmin, xmax)
ylab <- 'Normality test p-value'
ylim <- c(0, 0.1)
col1 <- "blue"
col2 <- "red"
cols <- c(col1, col2)

for (type in imageTypes){
    initializeImage("normtest", type)
    lines(xlim,c(0.05,0.05))
    plotShapiro <- function(ips, col){
        points(lens, ips, col=col, pch=20)
    }
    plotShapiro(ips,     col1)
    plotShapiro(ips_log, col2)
    legend("topleft",col=cols,lty=1,lwd=2,legend=c("unmodified","log10"),bg="white",cex=0.8)
}

# repeat input data to stdout
printCrosstab()

