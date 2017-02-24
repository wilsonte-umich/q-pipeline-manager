
# prepare the output directories
$Q_MOD_DIR      run echo \$Q_MOD_DIR
$LIB_DIR        $Q_MOD_DIR/bio/varBin/lib
$MAPS_DIR       $VAR_BIN_DIR/maps
$GROUP_DIR      $VAR_BIN_DIR/$GROUP_NAME
$BINS_DIR       $GROUP_DIR/bins
$HMM_DIR        $GROUP_DIR/HMM
invoke file/create.q $DIR $VAR_BIN_DIR $MAPS_DIR $GROUP_DIR $BINS_DIR $HMM_DIR

# set output file names
$REF_STAT_FILE  $BINS_DIR/$GROUP_NAME.ref.stats

# get the chromosome information
invoke bio/genome/chrom_info.q

preserve all
