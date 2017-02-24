#!/bin/bash

#q    require $SAMPLES $GROUP_NAME $REF_REGION $REF_PLOIDY
#q    require $LIB_DIR $MAPS_DIR $GROUP_DIR $BINS_DIR
#q    require $REF_STAT_FILE $GAP_FILE

#$    -N  vbSetCnt_$GROUP_NAME
#$    -wd $BINS_DIR
#$    -l  vf=4G

echo "optimizing $GROUP_NAME bin fragment count using $REF_REGION as reference"

perl $LIB_DIR/set_bin_count.pl $REF_REGION $REF_STAT_FILE
checkPipe

echo "done"
