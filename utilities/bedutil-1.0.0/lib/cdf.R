
#######################################################################
# 'cdf.R' does the work of command 'cdf'
#######################################################################

# collect and parse the input data
source(paste(Sys.getenv("SCRIPT_DIR"), "/common.R", sep=""))
loadCrosstab(TRUE)  # sorted by size
checkIterations()

# collect the simulated data for smallest, median and largest features
sm <- getFeature(1)            # sm = small
md <- getFeature(round(nF/2))  # md = median
bg <- getFeature(nF)           # bg = big

# common image properties
ylim <- c(0, 1)
col1 <- "blue"
col2 <- "red"
col3 <- "darkgreen"
cols <- c(col1, col2, col3)

#plot the data
rqSeries <- function(feat, col, discard){
    qnts <- sapply(feat$vs, getQuantile, feat)  # get quantiles for all iterations values
    points(1:feat$nv, qnts, col=col, pch=".")
}
cdfSeries <- function(feat, col, xs){
    lines(xs,  pnorm(xs, mean(feat$vs), sd(feat$vs)), col=col)
    points(feat$vs, (1:feat$nv)/feat$nv, col=col, pch=".")
}
createPlots <- function(suffix, plotSeries, xs){
    for (type in imageTypes){
        initializeImage(suffix, type)
        lines(xlim,c(0.5,0.5),lty=3)
        #lines(rep(xmax/2,2),ylim,lty=3)
        plotSeries(bg, col3, xs)  
        plotSeries(md, col2, xs)              
        plotSeries(sm, col1, xs)
        legend("topleft",col=cols,lty=1,lwd=2,legend=c(sm$l,md$l,bg$l),bg="white",cex=0.8)
        graphics.off()
    }
}
plotRQ <- function(){
    main <<- "Rank-Quantile Plot"
    ylab <<- 'Quantile' 
    xlab <<- 'Rank'
    xmin <<- 1
    xmax <<- max(sm$nv, md$nv, bg$nv)
    xlim <<- c(xmin, xmax)     
    createPlots("rank_quantile", rqSeries, NULL)
}
plotCDF <- function(suffix, label){
    main <<- "Cumulative Distribution Function"
    ylab <<- 'Cumulative Fraction' 
    xlab <<- label
    xmin <<- min(sm$vs, md$vs, bg$vs)
    xmax <<- max(sm$vs, md$vs, bg$vs)
    xlim <<- c(xmin, xmax)    
    xs <-  seq(xmin, xmax, length.out=1000)
    createPlots(suffix, cdfSeries, xs)    
}
plotRQ()
plotCDF("unmodified_cdf", env['SIM_LABEL'])
sm <- applyLog10(sm)
md <- applyLog10(md)
bg <- applyLog10(bg)
plotCDF("log10_cdf", paste("log10(", env['SIM_LABEL'], ")", sep=""))

# repeat input data to stdout
printCrosstab()

