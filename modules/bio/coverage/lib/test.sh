
function getTaskObject {
    local __VAR_NAME=$1
    local __INDEX=$(($TASK_ID + 1))
    local __OBJECT=${!__INDEX}
    eval $__VAR_NAME="'$__OBJECT'"
}

TASK_ID=2
getTaskObject VAR a b c
echo VAR=$VAR

