
#######################################################################
# 'by_score.R' does the work of command 'by_score'
#######################################################################

#----------------------------------------------------------------------
# collect data and set common variables
#----------------------------------------------------------------------

# define parameters
source(paste(Sys.getenv("SCRIPT_DIR"), "/common.R", sep=""))
sim_log <- asLogical(env['SIM_LOG'])  # 'sim' refers to the 2nd score assigned to features by 'crossing'
sim_log_prefix <- "" # ifelse(sim_log, "Log(", "")
sim_log_suffix <- "" # ifelse(sim_log, ")", "")
bed_log <- asLogical(env['BED_LOG'])  # 'bed' refers to the original BED feature score (field 5)
bed_log_prefix <- "" # ifelse(bed_log, "Log(", "")
bed_log_suffix <- "" # ifelse(bed_log, ")", "")
STRATUM_SIZE <- asNumeric(env['STRATUM_SIZE'])  # 'stratum' referes to rounded and limited BED scores
MIN_BED <- asNumeric(env['MIN_BED'])
MAX_BED <- asNumeric(env['MAX_BED'])
MIN_SIM <- asNumeric(env['MIN_SIM'])
MAX_SIM <- asNumeric(env['MAX_SIM'])
PEARSON <- asLogical(env['PEARSON'])
ADD_STRATA <- env['ADD_STRATA']
ADD_STRATA <- ifelse(ADD_STRATA == "", NA, ADD_STRATA)
COMPARE_STRATA <- env['COMPARE_STRATA']
COMPARE_STRATA <- ifelse(COMPARE_STRATA == "", NA, COMPARE_STRATA)
COMPARE_TEST <- env['COMPARE_TEST']
SUPPRESS_FIXED <- asLogical(env['SUPPRESS_FIXED'])
PERCENTILE <- asLogical(env['PERCENTILE'])
by_type <- ifelse(PERCENTILE, "by_percentile", "by_score")
if(PERCENTILE){
    bed_log <- FALSE
    bed_log_prefix <- "%ile "
    bed_log_suffix <- ""
}
STRATIN <- asNumeric(env['STRATIN'])
HEIGHT <- asNumeric(env['HEIGHT'])

# collect and parse the input data
loadCrosstab()

# common image properties
# TODO: could expose these and other plot properties as parameters
point_pch = 1
point_col = "grey"
agg_pch = 19
agg_col = "red"
qnt_col = "blue"
box_pch = 1
maxPoints <- 250  # boxplot is used if any stratum has >maxPoints
                  # all points omitted from raw score PDF (not jpg) if total _unique_ points >maxPoints
maxOut <- 100     # maximum number of _unique_ outlier points plotted per boxplot, i.e. per stratum
                  # will sample and show extreme values if there are more than this many unique outliers
stratin <- ifelse(!is.na(STRATIN), STRATIN, 0.25)    # x-axis width of one stratum, in inches (plot width is variable)
boxwex <- STRATUM_SIZE  # controls the width of boxes in boxplots
height <- ifelse(!is.na(HEIGHT), HEIGHT, 3)

#----------------------------------------------------------------------
# stratify BED features and correlate BED scores to simulation scores
#----------------------------------------------------------------------

stratifyBED(bed_log, sim_log)
corb <- ifelse(SUPPRESS_FIXED,    FALSE, nrow(cor$b$vs)>0)
corg <- ifelse(is.na(ADD_STRATA), FALSE, nrow(cor$g$vs)>0)

#----------------------------------------------------------------------
# aggregate simulation scores across BED score strata
#----------------------------------------------------------------------

agg <- aggregateStrata(sim_log)
aggb <- ifelse(nI==0 | SUPPRESS_FIXED,    FALSE, corb)
aggg <- ifelse(nI==0 | is.na(ADD_STRATA), FALSE, corg)

#----------------------------------------------------------------------
# perform pairwise group comparisons
#----------------------------------------------------------------------

compareGroups()

#----------------------------------------------------------------------
# initialize stats file
#----------------------------------------------------------------------

stats <- list(i=NULL,f=NULL)
STATS_PREFIX <- paste(env['IMAGE_PREFIX'], ".", by_type, sep="")
statsFile <- paste(STATS_PREFIX, ".stats.tsv", sep="")
cat("STATS_PREFIX\t", STATS_PREFIX, "\n", sep="", file=statsFile, append=FALSE)
cat("BED_LABEL\t", env['BED_LABEL'], "\n", sep="", file=statsFile, append=TRUE)
cat("SIM_LABEL\t", env['SIM_LABEL'], "\n", sep="", file=statsFile, append=TRUE)
cat("SIM_LOG\t", ifelse(sim_log, '+', ''), "\n", sep="", file=statsFile, append=TRUE)
cat("AGGREGATE_TYPE\t", agg_shrt, "\n", sep="", file=statsFile, append=TRUE)

#----------------------------------------------------------------------
# raw_score correlation plot
#----------------------------------------------------------------------

r <- getR("score", cor)

xlab <- paste(bed_log_prefix, env['BED_LABEL'], bed_log_suffix, sep="")  # raw_score always use log labeling when BED_LOG is true
xr <- setRawScoreXLim()

log_lab <- paste(sim_log_prefix, env['SIM_LABEL'], sim_log_suffix, sep="")
ylab <- log_lab
y <- setRawScoreYLim()

sprintfNumber <- function(n){
    if(n>=1e8) return( paste(round(n/1e9,1),"G",sep="") )
    if(n>=1e5) return( paste(round(n/1e6,1),"M",sep="") )
    if(n>=1e3) return( paste(round(n/1e3,1),"K",sep="") )
    n
}
limitRaw <- function(rw, min, max){
    rw <- cor$b$vs[[rw]]
    rw <- ifelse(rw > max, max, rw)
    rw <- ifelse(rw < min, min, rw)
    rw
}
raw <- data.frame(rwx=limitRaw("rwx",xmin,xmax), rwy=limitRaw("rwy",ymin,ymax))
raw <- aggregate(raw,list(rwx=raw$rwx,rwy=raw$rwy),length)
for (type in imageTypes){
    imageFile <- initializeImage(paste(by_type,".raw_score",sep=""), type)
    text(xr, y$top, labels=r)
    #if(type=="jpg" | nrow(raw)<=maxPoints){  # prevent slow-loading PDFs with too many data points
        points(raw[[1]],raw[[2]],pch=point_pch,col=point_col)
    #}
    points(cor$b$ag$st,cor$b$ag$x,pch=agg_pch,col=agg_col)
    lines(cor$b$ag$st,cor$b$ag$x,col=agg_col)
    nfs <- sapply(strata$b$nf, sprintfNumber)
    #text(strata$b$st, y$row2, labels=nfs)
    graphics.off()
    if(!is.na(PDF_CROP) & grepl(".pdf",imageFile)){
        system(paste("pdfcrop", imageFile, imageFile))
    }
}

#----------------------------------------------------------------------
# by_score plot
#----------------------------------------------------------------------

if(!PERCENTILE & asLogical(env['SUPPRESS_BED_LOG_LABEL'])){ xlab <- env['BED_LABEL'] }  # allow override of log() designation for users providing groups and group labels

xr <- setByScoreXLim()

allpts <- c(cor$b$vs$v, cor$g$vs$v, cor$b$ag$x, cor$g$ag$x)
ymin <- min(allpts, na.rm=TRUE)
ymax <- max(allpts, na.rm=TRUE)
y <- setYLim()

plotCor <- function(c, lines){  # the data points for actual value correlation plot
    if(strata$maxnf>maxPoints){
        d <- c$vs
        sts <- unlist(d$st)
        sts <- aggregate(sts,list(sts),length)[[1]]
        for(st in sts){
            myBoxplot(st,unlist(d[d$st==st,'v']),ocol=point_col)
        }
    } else {
        points(c$vs$st,c$vs$v,pch=point_pch,col=point_col)
    }
    points(c$ag$st,c$ag$x,pch=agg_pch,col=agg_col)
    if(lines){ lines(c$ag$st,c$ag$x,col=agg_col) }
}
divLine <- function(b, g){  # demarcation line between fixed-width and user-defined strata
    if(b & g){
        lines(rep(strata$maxb+STRATUM_SIZE/2,2),c(ymin,y$top),lty=3)
    }
}
xLabels <- function(){  # custom x-axis labels for strata
    axis(1, at=strata$sts, labels=strata$lbs)
}
nFeatures <- function(){  # a top label row indicating the number of features in each stratum
    nfs <- sapply(strata$nfs, sprintfNumber)
    text(strata$sts, y$row2, labels=nfs)
}
for (type in imageTypes){
    imageFile <- initializeImage(by_type, type, xaxt="n")
    text(xr, y$top, labels=r)  # report raw score correlation coefficient, even if bins not plotted
    if(corb){ plotCor(cor$b, TRUE)  }
    if(corg){
        plotCor(cor$g, FALSE)
        prevX <- NA
        prevY <- NA
        for(i in 1:length(grps)){
            st <- strata$maxb + i * STRATUM_SIZE
            ag <- unlist(cor$g$ag[cor$g$ag$st==st,"x"])
            if(length(ag) > 0){
                ag <- ag[1]
                if(i>1 & grpLineTo[i]){ lines(c(prevX,st), c(prevY,ag), col=agg_col) }
                prevX <- st
                prevY <- ag
            } else {
                prevX <- NA
                prevY <- NA
            }
        }
        add_pairwise("score", type, "xs")
    }
    xLabels()
    divLine(corb, corg)
    nFeatures()
    graphics.off()
    if(!is.na(PDF_CROP) & grepl(".pdf",imageFile)){
                system(paste("pdfcrop", imageFile, imageFile))
    }
}

#----------------------------------------------------------------------
# by_score.by_iteration plot
#----------------------------------------------------------------------

if(nI>0){

    ylab <- paste(agg_shrt, " ", log_lab, " ", sep="")  # other properties inherited from above
    ymin <- iagg_lim$min
    ymax <- iagg_lim$max
    y <- setYLim()

    for (type in imageTypes){
        imageFile <- initializeImage(paste(by_type,".by_iteration",sep=""), type, xaxt="n")
        add_bins(add_iteration, FALSE, agg_col, "i")
        add_groups(add_iteration, FALSE, agg_col, "i")
        xLabels()
        divLine(aggb, aggg)
        nFeatures()
        graphics.off()
        if(!is.na(PDF_CROP) & grepl(".pdf",imageFile)){
            system(paste("pdfcrop", imageFile, imageFile))
        }
    }

#----------------------------------------------------------------------
# by_score.by_feature plot
#----------------------------------------------------------------------

    if(agg$sgn){
        r <- getR("quantile", agg)

        ylab <- paste("Quantile(", log_lab, ")", sep="")
        ymin <- 0
        ymax <- 1
        y <- setYLim()

        for (type in imageTypes){
            imageFile <- initializeImage(paste(by_type,".by_feature",sep=""), type, xaxt="n")
            lines(xlim,rep(0,2))  # set boundary lines for quantile expectations
            lines(xlim,rep(0.5,2),lty=3)
            lines(xlim,rep(1,2))
            text(xr, y$top, labels=r)  # report quantile correlation coefficient, even if bins not plotted
            add_bins(add_feature, TRUE, qnt_col, "f")
            add_groups(add_feature, TRUE, qnt_col, "f")
            xLabels()
            divLine(aggb, aggg)
            add_pairwise("quantile", type, "qs")
            graphics.off()
            if(!is.na(PDF_CROP) & grepl(".pdf",imageFile)){
                system(paste("pdfcrop", imageFile, imageFile))
            }
        }
    }

}

#----------------------------------------------------------------------
# finish stats file
#----------------------------------------------------------------------

if(nI>0){  # simulated scores report actual aggregates and median quantiles, with p-values
    if(agg$sgn){ stats <- paste(stats$i, stats$f, sep=", ") } else { stats <- stats$i }
} else {   # non-simulated scores just report their aggregate value
    stats <- sprintf("%0.3g", strata$xs)
}
for(i in 1:length(strata$lbs)){
    cat(strata$lbs[[i]], "\t", stats[i], "\n", sep="", file=statsFile, append=TRUE)
}

#----------------------------------------------------------------------
# repeat input data to stdout
#----------------------------------------------------------------------

printCrosstab()


