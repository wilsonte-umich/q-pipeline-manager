
# get passed arguments
dataFile    <- Sys.getenv("dataFile")
eProbFile   <- Sys.getenv("eProbFile")
weight      <- as.numeric(Sys.getenv("weight"))

# declare variables
dCNs <- -2:2

# functions
getExpCov <- function(medCov, CN, dCN){
    weight * (medCov + medCov * dCN/CN)
}
eProbs <- function(v){
    medCov <- v[1]
    CN     <- v[2]
    rawCov <- v[3]
    binI   <- v[4]
    cns    <- dCNs + CN    
    if(CN == 0){
        exps <- ifelse(cns<0, 0, ifelse(cns==0, medCov,   medCov * (cns + 1)))
    } else {
        singleCopy <- getExpCov(medCov, CN, 1-CN)    
        noCopies   <- singleCopy * 0.05 + 1
        exps <- ifelse(cns<0, 0, ifelse(cns==0, noCopies, getExpCov(medCov, CN, cns-CN)))
    }
    paste(paste(dpois(round(rawCov,0), exps), collapse=","), binI, sep="\t")
}

# get data ($medCovCol, $copyNumCol, $rawCovCol, $binI)
data <- read.table(dataFile, header=FALSE, sep="\t")

# calculate bin-specific emission probabilities
write(apply(data, 1, eProbs), file=eProbFile, ncolumns=1, append=FALSE)
