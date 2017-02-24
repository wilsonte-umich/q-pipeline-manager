
#######################################################################
# 'aggregate.R' does the work of command 'aggregate'
#######################################################################

# collect and parse the input data
source(paste(Sys.getenv("SCRIPT_DIR"), "/common.R", sep=""))
loadCrosstab()
checkIterations()

# define variables
log <- asLogical(env['SIM_LOG'])
log_prefix <- ifelse(log, "Log(", "")
log_suffix <- ifelse(log, ")", "")

# histograms of row and columns aggregation results, all data
plotAggregate <- function(l, suffix, ref, p1lab, p2lab){
    if(!l$hist){
        cat("bedutil aggregate: data not valid for", suffix, "histogram plot\n", file=stderr())
        return()
    }

    xmin <<- min(l$av, l$mids, na.rm=TRUE)
    xmax <<- max(l$av, l$mids, na.rm=TRUE)
    xpad <- (xmax - xmin) * 0.1
    xlim <<- c(xmin-xpad, xmax+xpad)

    ylab <<- "Frequency"
    ymin <<- 0
    ymax <<- max(l$freqs)
    top <- ymax * 1.2
    ylim <<- c(ymin, top)
    yvert <- c(ymin, ymax * 1.15)

    barcol <- "blue"
    gauscol <- "red"
    actcol <- "black"

    cat("bedutil aggregate:", p1lab, "p-value:", sprintf("%0.2g",l$p1), "\n", file=stderr())
    cat("bedutil aggregate:", p2lab, "p-value:", sprintf("%0.2g",l$p2), "\n", file=stderr())
    p <- ifelse(!is.na(l$p2), l$p2, l$p1)
    p <- sprintf("%0.2g", p)
    plab <- ifelse(p==0, paste("p <", 1/nI), paste("p =", p))

    for (type in imageTypes){
        initializeImage(suffix, type)
        if(!is.na(ref)){ lines(rep(ref,2),yvert,lty=3) }
        points(l$mids,l$freqs,type='h',col=barcol)
        if(l$norm){
            x <- seq(from=xmin,to=xmax,length.out=1000)
            y <- l$a*(exp(-((x-l$m)^2)/(2*(l$s^2))))
            lines(x,y,col=gauscol)
        }
        if(!is.na(l$av)){
            lines(rep(l$av,2),yvert,col=actcol)
        }
        if(!is.na(p)){ text(ifelse(is.na(ref),l$av,ref), top, labels=plab ) }
        graphics.off()
    }
}

# run analyses and plot
d <- data
l <- aggregateFeatures(log=log)
suffix <- "by_feature"
main <- env['SIM_LABEL']
xlab <- "Quantile"
ref <- 0.5
p1lab <- "Sign"
p2lab <- "Wilcoxon"
discard <- plotAggregate(l, suffix, ref, p1lab, p2lab)

l <- aggregateIterations(log=log)
suffix <- "by_iteration"
main <- ""
xlab <- paste(agg_shrt, "(", log_prefix, env['SIM_LABEL'], log_suffix, ")", sep="")
ref <- NA
p1lab <- "Tail"
p2lab <- "Gaussian"
discard <- plotAggregate(l, suffix, ref, p1lab, p2lab)

# repeat input data to stdout
printCrosstab()

