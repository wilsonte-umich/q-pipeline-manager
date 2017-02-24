#!/bin/bash

#q    require $SAMPLES $GROUP_NAME $CHROMS $REF_PLOIDY $MAX_CN $IGNORE_CHROMS $MIN_CNV_BINS
#q    require $LIB_DIR $GROUP_DIR $BINS_DIR $HMM_DIR
#q    require $REF_STAT_FILE

#$    -N  vbAssign_$GROUP_NAME
#$    -wd $BINS_DIR
#$    -l  vf=2G

echo "assigning copy numbers for $GROUP_NAME"

HMM_FILE=$HMM_DIR/$GROUP_NAME.bins.hmm.bgz
PLOT_FILE=$HMM_DIR/$GROUP_NAME.plot.bgz
CNVS_FILE=$HMM_DIR/$GROUP_NAME.CNVs.txt

export TMP_DIR=/tmp/q_module_varBins
mkdir -p $TMP_DIR

perl $LIB_DIR/assign_copy_number.pl |
bgzip -c |
slurp -o $HMM_FILE
checkPipe

echo
slurp $HMM_FILE |
bgzip -cd |
perl $LIB_DIR/fill_plot_values.pl |
bgzip -c |
slurp -o $PLOT_FILE
checkPipe

tabix -p bed $PLOT_FILE
checkPipe

echo
slurp $PLOT_FILE |
bgzip -cd |
perl $LIB_DIR/find_CNVs.pl |
slurp -o $CNVS_FILE
checkPipe

echo "done"
