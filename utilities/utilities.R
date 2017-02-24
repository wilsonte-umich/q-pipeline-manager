
# This file contains functions useful to R worker scripts.
# Initialize them in your worker script by adding line:
#
# source(paste(Sys.getenv("Q_UTIL_DIR"), "/utilities.R", sep=""));

#####################################################################
# functions that interpret array job $TASK_ID
#--------------------------------------------------------------------
# usage:
#     getTaskID
#     getTaskObject LIST_NAME  [LIST_NAME refers to a space-delimited passed variable]
#--------------------------------------------------------------------
getTaskID <- function(){
    TASK_ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
    if(is.na(TASK_ID)){
        TASK_ID <- as.numeric(Sys.getenv("PBS_ARRAYID"))
    }
    return(TASK_ID)
}
#--------------------------------------------------------------------
getTaskObject <- function(LIST_NAME){  # e.g. getTaskObject("SAMPLES")
    LIST <- Sys.getenv(LIST_NAME)
    LIST <- strsplit(LIST, " ")
    LIST <- LIST[[1]]
    TASK_ID <- getTaskID()
    return(LIST[TASK_ID])
}
#--------------------------------------------------------------------
getEnvObject <- function(LIST_NAME, INDEX){  # e.g. getEnvObject("SAMPLES", 5)
    LIST <- Sys.getenv(LIST_NAME)
    LIST <- strsplit(LIST, " ")
    LIST <- LIST[[1]]
    return(LIST[INDEX])
}
#####################################################################

