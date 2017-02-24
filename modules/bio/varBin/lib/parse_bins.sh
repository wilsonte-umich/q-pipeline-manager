#!/bin/bash

#q    require $SAMPLES $GROUP_NAME $CHROMS $CHROM_SIZES $IGNORE_CHROMS
#q    require $LIB_DIR $MAPS_DIR $GROUP_DIR $BINS_DIR
#q    require $REF_STAT_FILE $GAP_FILE

#$    -N  vbParse_$GROUP_NAME
#$    -wd $BINS_DIR
#$    -l  vf=4G
#$    -t  1-$N_CHROMS

getTaskObject CHROM $CHROMS
getTaskObject CHROM_SIZE $CHROM_SIZES

echo "parsing variable width bins for $GROUP_NAME $CHROM"

BINS_FILE=$BINS_DIR/$GROUP_NAME.$CHROM.bins.bgz
STATS_FILE=$BINS_FILE.stats

export TMP_DIR=/tmp/q_module_varBins
mkdir -p $TMP_DIR

perl $LIB_DIR/parse_bins.pl $CHROM $CHROM_SIZE $STATS_FILE |
bgzip -c |
slurp -o $BINS_FILE
checkPipe

echo "done"
