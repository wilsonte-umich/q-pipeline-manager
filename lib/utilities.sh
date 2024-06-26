
# this file contains functions provided by q to parsed shell scripts:
#   checkPredecessors
#   getTaskID
#   checkTaskID
#   checkForData
#   checkPipe
#   snipStream
#   snipFile
#   getTaskFile

#####################################################################
# these functions are automatically called by every target script
#--------------------------------------------------------------------
function checkPredecessors {  # check whether predecessors timed out
    if [ "$Q_Q_TYPE" = "SGE" ]; then  # predecessor time-out checking only required for SGE
        if [ "$Q_PREDECESSORS" != "" ]; then 
            sleep 2  # give the prior job's log file a moment to finish being written    
            IFS=":" read -a JOB_IDS <<< "$Q_PREDECESSORS"  # read the array of predecessor jobIDs
            for JOB_ID in "${JOB_IDS[@]}"
            do
            	local FAILED="`grep -L 'q: exit_status:' $Q_LOG_DIR/*.o$JOB_ID* 2>/dev/null`"  # returns filenames that failed to record an exit status, i.e. timed out
            	if [ "$FAILED" != "" ]; then 
            	    echo "predecessor job $JOB_ID failed to report an exit status (it probably timed out)"
            	    exit 100
            	fi
            done  
        fi
    fi
}
function getTaskID {  # $TASK_ID is not set if this is not an array job
    if [[ "$PBS_ARRAYID" = "" && "$SGE_TASK_ID" = "" && "$SLURM_ARRAY_TASK_ID" = "" ]]; then
        TASK_NUMBER=1
        TASK_ID=""     
    elif [ ${PBS_ARRAYID+1} ]; then
        TASK_NUMBER=$PBS_ARRAYID
        TASK_ID="--task-id $PBS_ARRAYID"
    elif [ ${SGE_TASK_ID+1} ]; then
        TASK_NUMBER=$SGE_TASK_ID
        TASK_ID="--task-id $SGE_TASK_ID"
    elif [ ${SLURM_ARRAY_TASK_ID+1} ]; then
        TASK_NUMBER=$SLURM_ARRAY_TASK_ID
        TASK_ID="--task-id $SLURM_ARRAY_TASK_ID"
    fi
}
#####################################################################


#####################################################################
# functions which may be called by user to perform various process checks
#--------------------------------------------------------------------
# usage:
#     checkTaskID
#     checkForData "some command"
#     checkPipe
#     waitForFile $FILE [$TIME_OUT]
#--------------------------------------------------------------------
function checkTaskID {  # ensure that TASK_ID was set
    if [ -z $TASK_ID ]; then 
        echo "TASK_ID not set for array job"
        exit 100
    fi
}
function checkForData {  # ensure that a data stream will have at least one line of data
    local COMMAND="$1"
    if [ "$COMMAND" = "" ]; then 
        echo "checkForData error: system command not provided"
        exit 100
    fi
    local LINE_1="`$COMMAND | head -n1`"
    if [ "$LINE_1" = "" ]; then
        echo "no data; exiting quietly"
        exit 0
    fi
}
function checkPipe {  # ensure that all commands in a pipe had exit_status=0
    local PSS=${PIPESTATUS[*]}
    for PS in $PSS; do
       if [[ ( $PS > 0 ) ]]; then
           echo "pipe error: [$PSS]"
           exit 100
       fi;
    done   
}
function waitForFile {  # wait for a file to appear on the file system; default timeout=60 seconds
    local FILE="$1"
    local TIME_OUT="$2"
    if [ "$FILE" = "" ]; then 
        echo "waitForFile error: file not provided"
        exit 100
    fi
    if [ "$TIME_OUT" = "" ]; then 
        local TIME_OUT=60
    fi 
    local ELAPSED=0
    while [ ! -s $FILE ]
    do
        sleep 2;
        let "ELAPSED += 2"
        if [ "$ELAPSED" -gt "$TIME_OUT" ]; then
            echo "waitForFile error: $FILE not found after $TIME_OUT seconds"
            exit 100
        fi  
    done
}
function checkFileExists {  # verify non-empty file, or first of glob if called as checkFileExists $GLOB
    local FILE="$1"
    if [ "$FILE" = "" ]; then 
        echo "checkFileExists error: file not provided"
        exit 100
    fi
    if [ ! -s "$FILE" ]; then
        echo "file empty or not found on node "`hostname`
        echo $FILE
        exit 100
    fi
}
#####################################################################


#####################################################################
# functions which may be called by user to excerpt data output into log files
#--------------------------------------------------------------------
# usage:
#     snipStream "some command" [$MAX_LINES]
#     snipFile $FILE [$MAX_LINES]
#--------------------------------------------------------------------
function snipStream {  # excerpt an output stream, e.g. "gunzip -c $ZIP_FILE"
    local SNIP_STREAM="$1"
    if [ "$SNIP_STREAM" = "" ]; then 
        echo "snipStream error: no stream provided"
        exit 100
    fi
    __HEAD_STREAM="$SNIP_STREAM"
    __TAIL_STREAM="$SNIP_STREAM | tac"
    getSnipMaxLines $2
    echoSnipLines 
}
function snipFile {  # excerpt a text file (not a binary file)
    getSnipParameters $1 $2
    echoSnipLines    
}
#--------------------------------------------------------------------
# functions called by snipStream and snipFile, not for use in target scripts
#--------------------------------------------------------------------
function getSnipParameters {  # read the input arguments
    local SNIP_FILE="$1"
    if [ ! -f "$SNIP_FILE" ]; then 
        echo "snipFile error: file not found:"
        echo "$SNIP_FILE"
        exit 100
    fi
    echo
    echo $SNIP_FILE   
    __HEAD_STREAM="cat $SNIP_FILE"
    __TAIL_STREAM="tac $SNIP_FILE"
    getSnipMaxLines $2
}
function getSnipMaxLines {  # determine how many lines to print
    __MAX_LINES="$1"
    if [ "$__MAX_LINES" = "" ]; then 
        __MAX_LINES=10
    fi 
    __HALF_MAX_LINES=$(($__MAX_LINES/2))
}
function echoSnipLines {  # generate the output
    local N_LINES=`wc -l <($__TAIL_STREAM) | cut -d " " -f1`
    if [ "$N_LINES" -gt "$__MAX_LINES" ]; then
        head -n$__HALF_MAX_LINES <($__HEAD_STREAM)
        echo "..."
        head -n$__HALF_MAX_LINES <($__TAIL_STREAM) | tac
    else 
        $__HEAD_STREAM
    fi
    echo
}
#####################################################################


#####################################################################
# functions which may be called by user to interpret array $TASK_NUMBER
#--------------------------------------------------------------------
# usage:
#     getTaskFile $FILE_GLOB  [sets global variable $TASK_FILE]
#     getTaskObject $VAR_NAME object1 object2 ...  [sets global variable $VAR_NAME]
#--------------------------------------------------------------------
function getTaskFile {
    TASK_FILE=${!TASK_NUMBER};
}
function getTaskObject {  # example:  getTaskObject MY_CHUNK $CHUNKS
    local __VAR_NAME=$1
    local __INDEX=$(($TASK_NUMBER + 1))
    local __OBJECT=${!__INDEX}
    eval $__VAR_NAME="'$__OBJECT'"
}
#####################################################################


