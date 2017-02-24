
#######################################################################
# 'common.R' contains functions used in simulation analysis
#######################################################################

#----------------------------------------------------------------------
# parameter handling
#----------------------------------------------------------------------

# extract environment variables, includes all parameters
env <- Sys.getenv()
asLogical <- function(v){
    if(is.na(v)){ return(FALSE) }
    if(v=="" | v==0){ return(FALSE) }
    return(TRUE)
}
asNumeric <- function(v){
    if(is.na(v)){ return(NA) }
    if(v=="" | v=="FALSE"){ return(NA) }
    return(as.numeric(v))
}
replaceNA <- !asLogical(env['SUPPRESS_NULL'])  # TRUE if NA should be replaced with 0
PERCENTILE <- FALSE  # for scripts that don't use it
NO_JPG <- env['NO_JPG']
PDF_CROP <- env['PDF_CROP']
SUPPRESS_AGG <- env['SUPPRESS_AGG']
TIMES_100 <- env['TIMES_100']
if(is.na(NO_JPG)){
    imageTypes <- c("jpg", "pdf")
} else {
    imageTypes <- c("pdf")
}

#----------------------------------------------------------------------
# data stream handling
#----------------------------------------------------------------------

# define common data variables, set by loadCrosstab
input <- NULL  # unmodified input data
data <- NULL   # potentially modified data
d <- NULL      # working copy of filtered data
avc <- 7       # column with the actual value
ivcs <- NULL   # columns with the iteration values
nI <- 0        # number of iterations in input table
nF <- 0        # number of features in input table
lens <- NULL   # vector of feature lengths
lsum <- NA     # sum of lens
wgts <- NULL   # vector of length weights
scrs <- NULL   # vector of raw BED feature scores (might be logged, but never binned or limited)
               # or a vector of percentiles of scrs if PERCENTILE == TRUE

# collect and parse the input data (as formatted by crosstab)
featureStats <- function(){
    nF <<- nrow(data)
    lens <<- data[[3]] - data[[2]]
    lsum <<- sum(lens)
    wgts <<- lens / lsum
    scrs <<- data[[5]]
    if(PERCENTILE){
        l <- list(vs=scrs,nv=nF)
        scrs <<- sapply(scrs,getQuantile,l) * 100
    }
}
loadCrosstab <- function(orderBySize=FALSE){
    data <<- read.table("stdin", header=TRUE, sep="\t", stringsAsFactors=FALSE, na.strings=c("NA","NULL",""))
    nc <- ncol(data)
    if(!is.na(TIMES_100)){ data[,avc:nc] <<- data[,avc:nc] * 100 }
    input <<- data
    if(orderBySize){ data <<- data[order(data[[3]] - data[[2]]),] }
    if(nc > avc){
        ivcs <<- (avc+1):nc
        nI <<- length(ivcs)
    }
    featureStats()
    if (nF==0){
        cat("bedutil: fatal error: no input features\n", file=stderr())
        quit(save="no", status=1)
    }
}

# check for iterations for functions that need them
checkIterations <- function(){
    if(nI==0){
        cat("bedutil: no simulation iterations: exiting quietly\n", file=stderr())
        printCrosstab()
        quit(save="no")
    }
}

# print the unmodified input to STDOUT for streaming
printCrosstab <- function(){
    write.table(input, file="", quote=FALSE, sep="\t", na = "NULL", row.names=FALSE, col.names=TRUE)
}

#----------------------------------------------------------------------
# feature/iteration loading
#----------------------------------------------------------------------

# apply NA replacement with zeros xor NA purging from scores vector
maskNA <- function(vs, noPurge=FALSE){
    if(replaceNA){  # always replace NA with 0 if not SUPPRESS_NULL
        vs[is.na(vs)] <- 0
    } else if(!noPurge) {  # allow caller to override SUPPRESS_NULL, required for actual value
        vs <- vs[!is.na(vs)]
    }
    vs
}

# apply log10 to actual and iteration score values, by feature row or iteration column
applyLog10 <- function(l, log=TRUE, type="feat"){
    if(log){
        l$av <- ifelse(l$av>0, log10(l$av), NA)  # actual can be NA
        if(type=="feat"){
            l$vs <- log10(l$vs[l$vs>0])          # feature iterations cannot be NA
        } else {
            l$vs[l$vs<=0] <- NA                  # iteration vectors keep all features, including NA
            l$vs <- log10(l$vs)
        }
        l$nv <- length(l$vs[!is.na(l$vs)])
    }
    l
}

# collect a specific iteration from data, by column index
getIteration <- function(col, log=FALSE){
    if(is.null(d)){ d <- data }
    vs <- maskNA(unlist(d[,col]), noPurge=TRUE)  # iteration vectors keep all features, including NA
    nv <- length(vs[!is.na(vs)])
    l <- list(av=NA,vs=vs,nv=nv,wgts=d[[3]]-d[[2]],wgts5=d[[5]])
    l <- applyLog10(l, log, "iter")
    l
}

# collect a specific feature from data, by row index
getFeature <- function(row, log=FALSE){
    if(is.null(d)){ d <- data }
    av <- maskNA(d[row,avc], noPurge=TRUE)  # actual value can be NA
    vs <- NA
    nv <- 0
    if(nI > 0){
        vs <- maskNA(unlist(d[row,ivcs]))  # feature iterations cannot be NA
        vs <- sort(vs)  # always sort a feature's iterations
        nv <- length(vs)
    }
    l <- list(av=av,vs=vs,nv=nv,l=lens[row],w=wgts[row])
    l <- applyLog10(l, log, "feat")
    l
}

#----------------------------------------------------------------------
# score aggregation
#----------------------------------------------------------------------

# create a histogram of aggregated values
scoreHist <- function(l){
    l$hist <- FALSE
    l$norm <- FALSE
    l$vs <- l$vs[!is.na(l$vs)]  # drop NA from aggregate lists
    l$nv <- length(l$vs)
    if(l$nv == 0){ return(l) }
    h <- hist(l$vs, as.numeric(env['N_BINS']), plot=FALSE)
    l$hist <- TRUE
    l$freqs <- h$counts/l$nv
    l$mids <- h$mids
    vs <- l$vs
    if(l$nv>5000){ vs <- sample(l$vs, 5000) }
    vs <- sort(vs)
    sh <- ifelse(l$nv>=3 & vs[1]!=vs[length(vs)], shapiro.test(vs), 0)
    if(sh>=0.05){  # only fit gaussian if data are ~normal
        l$a <- max(l$freqs)
        l$m <- sum(l$freqs * l$mids)
        freqSum <- cumsum(l$freqs)
        l$s <- min(l$mids[freqSum>=0.841]) - l$m
        guassian <- freqs ~ a*(exp(-((mids-m)^2)/(2*(s^2))))
        curveFit <- try(nls(guassian,data=data.frame(freqs=l$freqs,mids=l$mids),start=list(a=l$a,m=l$m,s=l$s)), silent=TRUE)
        curveFitOK <- !inherits(curveFit, "try-error")
        if(curveFitOK){
            l$norm <- TRUE
            coef <- coef(curveFit)
            l$a <- coef[['a']]
            l$m <- coef[['m']]
            l$s <- coef[['s']]
        }
    }
    l
}

# aggregate iterations
agg_cmd <- NULL
agg_shrt <- NULL
weighted <- NULL
weighted5 <- NULL
bool <- function(vs, na.rm=TRUE){ # non-standard aggregate function
    if(na.rm){ vs <- vs[!is.na(vs)] }
    length(vs[vs!=0&vs!=""])
}
set_agg <- function(){  # interpret the requested aggregate type
    agg_cmd <<- match.fun(switch(env['AGGREGATE_TYPE'],
        'AVERAGE'   =  'mean',
        'WEIGHTED'  =  'weighted.mean',
        'WEIGHTED5' =  'weighted.mean',
        'SUM'       =  'sum',
        'COUNT'     =  'bool',
        'MEDIAN'    =  'median',
        'MIN'       =  'min',
        'MAX'       =  'max'
    ))
    agg_shrt <<- switch(env['AGGREGATE_TYPE'],
        'AVERAGE'   =  'Avg.',
        'WEIGHTED'  =  'Wt. Avg.',
        'WEIGHTED5' =  'Wt. Avg.',
        'SUM'       =  'Sum',
        'COUNT'     =  'Count',
        'MEDIAN'    =  'Median',
        'MIN'       =  'Min',
        'MAX'       =  'Max'
    )
    if(!is.na(SUPPRESS_AGG)){agg_shrt <<- ""}
    weighted  <<- env['AGGREGATE_TYPE'] == 'WEIGHTED'
    weighted5 <<- env['AGGREGATE_TYPE'] == 'WEIGHTED5'
}
aggregateIteration <- function(col, log){  # returns an iteration aggregate
    l <- getIteration(col, log)
    if(weighted){
        weighted.mean(l$vs, l$wgts, na.rm=TRUE)
    } else if(weighted5){
        weighted.mean(l$vs, l$wgts5, na.rm=TRUE)
    } else {
        agg_cmd(l$vs, na.rm=TRUE)
    }
}
aggregateIterations <- function(log=FALSE){  # compare actual aggregate score against histogram of iteration aggregates
    set_agg()
    l <- list(av=aggregateIteration(avc, log))
    l$vs <- sapply(ivcs, aggregateIteration, log)
    l <- scoreHist(l)
    l$p1 <- NA
    l$p2 <- NA
    if(is.na(l$av) | l$nv == 0){ return(l) }
    l$p1 <- min(length(l$vs[l$vs<=l$av]), length(l$vs[l$vs>=l$av])) / l$nv  # one-sided distribution tail
    if(!l$norm){ return(l) }
    l$p2 <- pnorm(l$av,mean=l$m,sd=l$s,lower.tail=(l$av<=l$m)) * 2  # two-sided normal p
    l
}

# aggregate features
sgn <- TRUE    # set to NA if sign test is invalid for any feature
wcx <- TRUE    # set to NA if Wilcoxon test is invalid for any feature
getQuantile <- function(v, l){  # returns the quantile of a provided score among all iteration scores
    (length(l$vs[l$vs<v]) + length(l$vs[l$vs==v]) / 2) / l$nv
}
aggregateFeature <- function(row, log){  # returns the iteration quantile of a feature's actual score
    if(is.na(sgn)){ return(NA) }  # no point in continuing to analyze features, quantiles meaningless for this score
    l <- getFeature(row, log)
    if(is.na(l$av) | l$nv<3){ return(NA) }
    med <- median(l$vs)
    min <- min(l$vs)
    max <- max(l$vs)
    nmin <- length(l$vs[l$vs==min])
    nmax <- length(l$vs[l$vs==max])
    threshold <- l$nv/20
    if(med==min | med==max){  # sign test (and thus Wilcoxon) requires that median is uniquely defined as different than extremes
        sgn <<- NA
        wcx <<- NA
        return(NA)
    } else if(nmin>threshold | nmax>threshold){  # reject bounded score types for Wilcoxon
        wcx <<- NA                               # TODO: elaborate the Wilcoxon criteria?
    }
    getQuantile(l$av, l)
}
aggregateFeatures <- function(log=FALSE){  # compare the distribution of feature score quantiles to expected median of 0.5
    l <- list(av=NA)
    l$vs <- sapply(1:nrow(d), aggregateFeature, log)
    l$allvs <- l$vs  # will contain NA quantiles, for correlation to scrs
    l <- scoreHist(l)
    l$p1 <- NA
    l$p2 <- NA
    if(!is.na(sgn)){  # apply the sign test
        sgns <- sign(l$vs - 0.5)
        npos <- length(sgns[sgns>0])
        nv <- length(sgns[sgns!=0])
        l$p1 <- binom.test(npos, nv, p=0.5)$p.value
    }
    if(!is.na(wcx)){  # apply the Wilcoxon test
        l$p2 <- wilcox.test(l$vs, mu=0.5)$p.value
    }
    l
}

#----------------------------------------------------------------------
# BED feature score stratification
#----------------------------------------------------------------------

strata <- list(mins=NA,maxb=0,maxs=NA,sts=NULL,lbs=NULL,nfs=NULL,maxnf=NULL,xs=NULL,qs=NULL)  # strata define the BED-score groupings on the x-axis of by_score plots
grps <- NULL   # a vector of user-requested strata
grpLineTo <- NULL # a vector of booleans whether to draw a line to the group aggregate point
cor <- list()  # holds correlations of BED scores to actual second/simulation scores
iagg_lim <- list(min=NA,max=NA)
logDataFrame <- function(d, log, col){  # apply log10 to a specific column of a data frame, otherwise apply NA handling
    if(log){
        d <- d[d[[col]]>0,]
        d[,col] <- log10(d[,col])
        d <- d[!is.na(d[[col]]),]
    } else if(replaceNA){
        d[is.na(d[[col]]),col] <- 0
    } else {
        d <- d[!is.na(d[[col]]),]  # NA scores are useless features for by_score
    }
    d
}
limitBED <- function(d, col){  # push BED scores into strata
    bin_fxn <- ifelse(PERCENTILE, "floor", "round")
    bin_fxn <- match.fun(bin_fxn)
    d[,col] <- bin_fxn(d[,col]/STRATUM_SIZE)*STRATUM_SIZE  # bin the BED scores
    if(!is.na(MIN_BED)){ d[,col] <- sapply(d[[col]], max, MIN_BED) } # apply strata limits
    if(!is.na(MAX_BED)){ d[,col] <- sapply(d[[col]], min, MAX_BED) }
    d
}
limitSIM <- function(d, col){  # apply plot limits to simulation scores
    if(!is.na(MIN_SIM)){ d[,col] <- sapply(d[[col]], max, MIN_SIM) }
    if(!is.na(MAX_SIM)){ d[,col] <- sapply(d[[col]], min, MAX_SIM) }
    d
}
between <- function(x,v){  # operator for semi-inclusive between comparisons in user-defined strata
    x>=v[1] & x<v[2]
}
parseGrp <- function(i, scrs, df=NULL){  # assemble user-defined stratum
    g <- list()
    g$st <- strata$maxb + i * STRATUM_SIZE
    grp <- strsplit(grps[i],',')[[1]]
    lb <- grp[1]
    grpLineTo <<- c(grpLineTo, length(grep('^!',lb))==0)
    g$lb <- sub('^!','',lb)
    op <- match.fun(grp[2])
    comps <- as.numeric(grp[3:length(grp)])
    if(is.null(df)){ df <- data }
    is <- op(scrs,comps)
    d <<- df[is,]  # hold stratum data in d global working variable
    g$nf <- nrow(d)
    g
}
stratifyBED <- function(bed_log=FALSE, sim_log=FALSE){  # generate strata for by_score correlation plot of bed score vs. actual simulation score

    # apply log to BED feature score
    data <<- logDataFrame(data, bed_log, 5)
    featureStats() # reset lens, wgts and scrs based on logged data

    # prepare a new tmp data frame with columns of logged and logged and limited actual data
    df <- data.frame(lbs=data[[5]], lss=data[[avc]], wgts=wgts)
    df <- logDataFrame(df, sim_log, 'lss')  # might remove rows
    if(PERCENTILE){
        l <- list(vs=df$lbs,nv=nrow(df))
        df$lbs <- sapply(df$lbs,getQuantile,l) * 100
    }
    df$llbs <- df$lbs
    df$llss <- df$lss
    df <- limitBED(df, 'llbs')
    df <- limitSIM(df, 'llss')
    set_agg()

    # calculate correlation coefficients and trendline on logged but not rounded or limited data
    cor$rp <<- round( cor(df$lbs, df$lss, method="pearson")  ,3)
    cor$rs <<- round( cor(df$lbs, df$lss, method="spearman") ,3)

    # fixed-width strata ('bins')
    # always do this as it establishes bin strata required for quantile correlation analysis
        cor$b$vs <<- data.frame(rwx=df$lbs,rwy=df$lss,st=df$llbs,v=df$llss)
        bs <- list(st=df$llbs)
        if(weighted){  # aggregate the actual features by bin, prior to simulation score limit
            wsums <- aggregate(df$wgts, bs, sum, na.rm=TRUE)$x
            cor$b$ag <<- aggregate(df$lss*df$wgts, bs, sum, na.rm=TRUE)
            cor$b$ag$x <<- cor$b$ag$x / wsums
        } else if(weighted5){
            wsums <- aggregate(df$lbs, bs, sum, na.rm=TRUE)$x
            cor$b$ag <<- aggregate(df$lss*df$lbs, bs, sum, na.rm=TRUE)
            cor$b$ag$x <<- cor$b$ag$x / wsums
        } else {
            cor$b$ag <<- aggregate(df$lss, bs, agg_cmd, na.rm=TRUE)
        }
        strata$b$st <<- unlist(cor$b$ag$st)
        strata$b$nf <<- unlist(aggregate(df$lss, bs, length)$x)
        strata$b$lb <<- strata$b$st
    if(!SUPPRESS_FIXED){  # but only record the bins in plotting strata when told to do so
        strata$mins <<- min(strata$b$st, na.rm=TRUE)
        strata$maxb <<- max(strata$b$st, na.rm=TRUE)
        strata$maxs <<- strata$maxb
        strata$sts <<- strata$b$st
        strata$lbs <<- strata$b$lb
        strata$lbs <<- ifelse(!is.na(MIN_BED) & strata$lbs==MIN_BED, paste("<",strata$lbs+STRATUM_SIZE,sep=""), strata$lbs)
        strata$lbs <<- ifelse(!is.na(MAX_BED) & strata$lbs==MAX_BED, paste(">",strata$lbs-STRATUM_SIZE,sep=""), strata$lbs)
        strata$nfs <<- strata$b$nf
        strata$maxnf <<- max(strata$nfs, na.rm=TRUE)
        strata$xs <<- unlist(cor$b$ag$x)
    }

    if(!is.na(ADD_STRATA)){  # user-defined strata ('groups')
        grps <<- strsplit(ADD_STRATA,';')[[1]]
        cor$g$vs <<- data.frame(st=NULL,v=NULL)
        cor$g$ag <<- data.frame(st=NULL,x=NULL)
        for(i in 1:length(grps)){
            g <- parseGrp(i, df$lbs, df)
            if(g$nf>0){
                v <- ifelse(weighted, weighted.mean(d$lss, d$wgts, na.rm=TRUE), ifelse(weighted5, weighted.mean(d$lss, d$lbs, na.rm=TRUE), agg_cmd(d$lss, na.rm=TRUE)))
                cor$g$vs <<- rbind(cor$g$vs, data.frame(st=g$st,v=d$llss))
                cor$g$ag <<- rbind(cor$g$ag, data.frame(st=g$st,x=v))
                strata$g$st <<- c(strata$g$st, g$st)
                strata$g$nf <<- c(strata$g$nf, g$nf)
                strata$g$lb <<- c(strata$g$lb, g$lb)
            } else {
                cat("bedutil by_score: stratifyBED: group", g$lb, "did not contain any features\n", file=stderr())
            }
        }
        strata$mins <<- min(strata$mins, strata$g$st, na.rm=TRUE)
        strata$maxs <<- max(strata$maxs, strata$g$st, na.rm=TRUE)
        strata$sts <<- c(strata$sts, strata$g$st)
        strata$lbs <<- c(strata$lbs, strata$g$lb)
        strata$nfs <<- c(strata$nfs, strata$g$nf)
        strata$maxnf <<- max(strata$maxnf, strata$nfs, na.rm=TRUE)
        strata$xs <<- c(strata$xs, unlist(cor$g$ag$x))
    }

}

aggregateStratum <- function(log){  # apply feature and iteration aggregation to a BED score stratum
    l <- list(test="TEST")
    sgn <<- TRUE  # process sign and Wilcoxon test permissions indepently for each stratum
    wcx <<- TRUE
    l$i <- aggregateIterations(log)
    l$f <- aggregateFeatures(log)
    l$f$sgn <- sgn
    l$f$wcx <- wcx
    l
}
aggregateStrata <- function(log=FALSE){  # generate strata for by_iteration and by_feature enrichment plots

    # round and limit the BED feature scores
    if(PERCENTILE){ data[,5] <<- scrs }
    data <<- limitBED(data, 5)

    agg <- list(sgn=FALSE,b=list(),g=list())  # for unknown reasons, this fails when defined globally!
                                              # must define locally and pass list back to global environment

    if(nI>0){  # false if not a simulation set, i.e. no data to process

        # fixed-width strata ('bins')
        # always do this so that quantile correlation is established
            cr <- data.frame(st=NULL,vs=NULL)
            for (b in strata$b$st){
                is <- data[[5]]==b
                d <<- data[is,]
                bc <- as.character(b)
                agg$b[[bc]] <- aggregateStratum(log)
                agg$b[[bc]]$nf <- strata$b$nf[strata$b$st==b]
                agg$b[[bc]]$st <- b
                agg$b[[bc]]$lb <- b
                agg$sgn <- agg$sgn | !is.na(sgn)
                cr <- rbind(cr, data.frame(st=scrs[is], vs=agg$b[[bc]]$f$allvs))  # correlation uses unbinned BED scores and feature quantiles
                if(!SUPPRESS_FIXED){
                    iagg_lim$min <<- min(iagg_lim$min, agg$b[[bc]]$i$av, agg$b[[bc]]$i$vs, na.rm=TRUE)
                    iagg_lim$max <<- max(iagg_lim$max, agg$b[[bc]]$i$av, agg$b[[bc]]$i$vs, na.rm=TRUE)
                }
            }
            cr <- cr[!is.na(cr$st)&!is.na(cr$vs),]
            agg$rp <- round( cor(cr$st, cr$vs, method="pearson")  ,3)
            agg$rs <- round( cor(cr$st, cr$vs, method="spearman") ,3)

        if(!is.na(ADD_STRATA)){  # user-defined strata ('groups')
            for(i in 1:length(grps)){
                g <- parseGrp(i, scrs)
                ic <- as.character(i)
                if(g$nf>0){
                    agg$exists[[ic]] <- TRUE
                    agg$g[[ic]] <- aggregateStratum(log)
                    agg$g[[ic]]$nf <- strata$g$nf[i]
                    agg$g[[ic]]$st <- g$st
                    agg$g[[ic]]$lb <- g$lb
                    agg$sgn <- agg$sgn | !is.na(sgn)
                    iagg_lim$min <<- min(iagg_lim$min, agg$g[[ic]]$i$av, agg$g[[ic]]$i$vs, na.rm=TRUE)
                    iagg_lim$max <<- max(iagg_lim$max, agg$g[[ic]]$i$av, agg$g[[ic]]$i$vs, na.rm=TRUE)
                } else {
                    cat("bedutil by_score: aggregateStrata: group", g$lb, "did not contain any features\n", file=stderr())
                    agg$exists[[ic]] <- FALSE
                }
            }
        }

    }

    agg
}

prs <- list()
getPairP <- function(vs1,vs2,test){
    vs1 <- vs1[!is.na(vs1)]
    vs2 <- vs2[!is.na(vs2)]
    p <- NA
    if(length(vs1)<2 | length(vs2)<2){
        # do nothing
    } else if(test=="t.test"){
        p <- t.test(vs1,vs2)$p.value
    } else if (test=="wilcox.test"){
        p <- wilcox.test(vs1,vs2,exact=FALSE)$p.value
    } else if (test=="fisher.test"){
        n1 <- length(vs1[vs1==0])
        p1 <- length(vs1[vs1==1])
        n2 <- length(vs2[vs2==0])
        p2 <- length(vs2[vs2==1])
        m <- matrix(c(n1,p1,n2,p2),ncol=2)
        p <- fisher.test(m)$p.value
    } else {
        cat("bedutil: fatal error: unknown COMPARE_TEST:", test, "\n", file=stderr())
        quit(save="no", status=1)
    }
    return(p)
}
compareGroups <- function(){  # perform pairwise group comparisons
    if(is.na(COMPARE_STRATA)){ return() }
    prs_ <- strsplit(COMPARE_STRATA,';')[[1]]
    for(i in 1:length(prs_)){
        ic <- as.character(i)
        pr <- strsplit(prs_[i],',')[[1]]
        lb1 <- pr[1]
        lb2 <- pr[2]
        i1 <- match(lb1,strata$g$lb)
        i2 <- match(lb2,strata$g$lb)
        st1 <- strata$g$st[i1]
        st2 <- strata$g$st[i2]
        vs1 <- unlist(cor$g$vs[cor$g$vs$st==st1,"v"])
        vs2 <- unlist(cor$g$vs[cor$g$vs$st==st2,"v"])
        ic1 <- as.character(i1)
        ic2 <- as.character(i2)
        qs1 <- agg$g[[ic1]]$f$vs
        qs2 <- agg$g[[ic2]]$f$vs
        vp <- getPairP(vs1,vs2,COMPARE_TEST)
        qp <- getPairP(qs1,qs2,"wilcox.test")
        prs[[ic]]$pr <<- prs_[i]
        prs[[ic]]$st1 <<- st1
        prs[[ic]]$st2 <<- st2
        prs[[ic]]$ps <<- list(score=vp,quantile=qp)
    }
}

#----------------------------------------------------------------------
# image/plot handling
#----------------------------------------------------------------------

# common image properties
width <- 3.5
height <- 3.167
units <- 'in' #w and h in inches
pointsize <- 8
res <- 600 #dpi
pch <- 20
cex <- 0.3
main <- ""
xmin <- NULL  # limits and labels must be set by calling script prior to initialization
xmax <- NULL
xlim <- NULL
xlab <- ""
ymin <- NULL
ymax <- NULL
ylim <- NULL
ylab <- ""

# initialize plot image
initializeImage <- function(suffix, type, xaxt=NULL){
    imageFile <- paste(env['IMAGE_PREFIX'], ".", suffix, ".", type, sep="")
    imageFile <- gsub(" ", "_", imageFile)
    if(type == "jpg"){
        jpeg(file=imageFile,width=width,height=height,unit=units,pointsize=pointsize,res=res)
    } else {
        pdf(file=imageFile,width=width,height=height,pointsize=pointsize)
    }
    par(mgp=c(2.5,1,0))
    plot(xlim,ylim,type="n",main=main,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,xaxt=xaxt)
    imageFile
}

# parse p-values and correlation coefficients for plot labeling
getSig <- function(p1, p2){
    p <- ifelse(is.na(p2), p1, p2)
    pp <- ""
    if(is.na(p)){ pp <- "" }
    else if(p<0.0001){ pp <- "***" }
    else if(p<0.001){ pp <- "**" }
    else if(p<0.01){ pp <- "*" }
    return(pp)
}
getR <- function(type, l){
    cat("bedutil", by_type, type, "Pearson:", l$rp, "\n", file=stderr())
    cat("bedutil", by_type, type, "Spearman:", l$rs, "\n", file=stderr())
    r <- ifelse(PEARSON, l$rp, l$rs)
    cat(paste(type, "r", sep=" "), "\t", r, "\n", sep="", file=statsFile, append=TRUE)
    paste("r = ", r, sep="")
}

# axis functions
setRawScoreXLim <- function(){
    xmin <<- max(min(cor$b$vs$rwx), MIN_BED, na.rm=TRUE)
    xmax <<- min(max(cor$b$vs$rwx), MAX_BED, na.rm=TRUE)
    xr <- xmin + (xmax - xmin)/2
    width <<- 3.5
    xlim <<- c(xmin, xmax)
    xr
}
setRawScoreYLim <- function(){
    ymin <<- max(min(cor$b$vs$rwy), MIN_SIM, na.rm=TRUE)
    ymax <<- min(max(cor$b$vs$rwy), MAX_SIM, na.rm=TRUE)
    setYLim()
}
setByScoreXLim <- function(){
    nst <- ((strata$maxs - strata$mins) / STRATUM_SIZE) + 1
    xrms <- ifelse(SUPPRESS_FIXED, strata$maxs, strata$maxb)
    xr <- strata$mins + (xrms - strata$mins)/2
    xmin <<- strata$mins
    xmax <<- strata$maxs
    xin <- nst * stratin  # includes half stratum padding on each end
    width <<- xin + par("mai")[2] + par("mai")[4]
    xlim <<- c(xmin - STRATUM_SIZE/2, xmax + STRATUM_SIZE/2)
    xr
}
setYLim <- function(){
    y <- list()
    y$wid <- (ymax - ymin)
    if(y$wid==0){
        ymin <<- ymin - 1
        ymax <<- ymax + 1
        y$wid <- (ymax - ymin)
    }
    y$top <- ymin + y$wid * 1.2
    ylim <<- c(ymin, y$top)
    y$row2 <- ymin + y$wid * 1.085
    y
}

# functions that plot by_score elements
myBoxplot <- function(x,ys,ocol="black"){  # boxplot with outlier point number control
    bp <- boxplot(ys,add=TRUE,at=x,outline=FALSE,boxwex=boxwex,yaxt="n")
    nout <- length(bp$out)
    if(nout>0){
        out <- aggregate(bp$out,list(bpo=bp$out),length)[[1]]
        nout <- length(out)
        if(nout>maxOut){
            omin <- min(out, na.rm=TRUE)  # be sure to plot the extremes
            omax <- max(out, na.rm=TRUE)
            out <- c(omin,omax,sample(out,maxOut))
        }
        points(rep(x,length(out)),out,pch=box_pch,col=ocol)
    }
}
add_iteration <- function(a){  # aggregated by iteration, i.e. crosstab column
    myBoxplot(unlist(a$st),unlist(a$i$vs))
    text(a$st, y$top, labels=getSig(a$i$p1, a$i$p2))
}
add_feature <- function(a){  # aggregated by feature, i.e. crosstab row
    if(!is.na(a$f$sgn)){
        if(strata$maxnf>maxPoints){
            myBoxplot(unlist(a$st),unlist(a$f$vs),ocol=point_col)
        } else {
            points(rep(a$st,a$nf),a$f$vs,pch=point_pch,col=point_col)
        }
        text(a$st, y$row2, labels=getSig(a$f$p1, a$f$p2))
    }
}
getAV <- function(a, isFeat, type){  # get aggregated value and set related entry for stats file
    av <- ifelse(isFeat, median(a[[type]]$vs,na.rm=TRUE), a[[type]]$av)  # always aggregate quantiles as medians
    if(isFeat){ strata$qs <<- c(strata$qs, av) }
    p <- ifelse(is.na(a[[type]]$p2), a[[type]]$p1, a[[type]]$p2)
    stats[[type]] <<- c(stats[[type]], paste(sprintf("%0.3g",av), " (", sprintf("%0.3g",p), ")", sep=""))
    av
}
add_bins <- function(add_stratum, isFeat, agg_col, type){  # plot fixed-width strata ('bins')
    if(aggb){
        avs <- NULL
        for (b in strata$b$st){
            bc <- as.character(b)
            a <- agg$b[[bc]]
            add_stratum(a)
            av <- getAV(a, isFeat, type)
            avs <- c(avs, av)  # always aggregate quantiles as medians
        }
        points(strata$b$st,avs,pch=agg_pch,col=agg_col)
        lines(strata$b$st,avs,col=agg_col)  # only fixed-width, i.e. unit-spaced, strata aggregates are connected by lines
    }
}
add_groups <- function(add_stratum, isFeat, agg_col, type){  # plot user-defined strata ('groups')
    if(aggg){
        prevX <- NA
        prevY <- NA
        for(i in 1:length(grps)){
            ic <- as.character(i)
            if(agg$exists[[ic]]){
                a <- agg$g[[ic]]
                add_stratum(a)
                av <- getAV(a, isFeat, type)
                points(a$st, av, pch=agg_pch, col=agg_col)
                if(i>1 & grpLineTo[i]){ lines(c(prevX,a$st), c(prevY,av), col=agg_col) }
                prevX <- a$st
                prevY <- av
            } else {
                prevX <- NA
                prevY <- NA
            }
        }
    }
}

# functions that plot pairwise score comparisons
add_pairwise <- function(type, imageType, aggs){
    if(is.na(COMPARE_STRATA)){ return() }
    yinc <- y$wid/15  # determine the spacing of pairwise bars
    minagg <- min(strata[[aggs]], na.rm=TRUE)  # determine whether there is more room on top or bottom for bars
    maxagg <- max(strata[[aggs]], na.rm=TRUE)
    topwid <- ymax - maxagg
    botwid <- minagg - ymin
    ptop <- ifelse(topwid>=botwid, ymax, minagg)  # position the bars at plot top or just below the lowest aggregate value
    j <- 0
    for(i in 1:length(prs)){
        ic <- as.character(i)
        st1 <- prs[[ic]]$st1
        st2 <- prs[[ic]]$st2
        p <- prs[[ic]]$ps[[type]]
        if(imageType=="jpg"){  # always report the pairwise p values
            stat <- paste(prs[[ic]]$pr, type, "p-value")
            pf <- sprintf("%0.3g",p)
            cat("bedutil by_score", stat, ":", pf, "\n", file=stderr())
            cat(stat, "\t", pf, "\n", sep="", file=statsFile, append=TRUE)
        }
        if(!is.na(p) & p<0.01){  # but only plot if significant
            py <- ptop - yinc*(1+j)
            lines(c(st1,st2),c(py,py),lwd=0.5)
            text(st1+STRATUM_SIZE/2, py-yinc*0.4, labels=getSig(p, NA))
            j <- j + 1
        }
    }
}



