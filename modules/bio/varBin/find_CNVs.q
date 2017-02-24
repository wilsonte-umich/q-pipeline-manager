
# set the common variables and paths
invoke lib/common.q
$IGNORE_CHROMS  $IGNORE_CHROMS ? $IGNORE_CHROMS : NO_DATA
$IGNORE_CHROMS  $IGNORE_CHROMS, # helps with pattern matching

# determine and set the optimal value of $VAR_BIN_COUNT
qsub lib/set_bin_count.sh

# parse to equally-weighted variable-width bins across all input samples
qsub lib/parse_bins.sh

# determine the modal copy number across all samples by Hidden Markov and find runs
$MAX_CN     $MAX_CN ? $MAX_CN : 4
qsub lib/assign_copy_number.sh








exit


# calls CNVs for each sample by looking for bin runs different from samples median
$CALL_DIR   $CALL_DIR/$GROUP
invoke file/create.q $DIR $CALL_DIR
invoke $MASTERS_DIR/coverage/slaves/callCNVs.q $SAMPLE $SAMPLES

# assemble final graphical output
$N_SAMPLES run echo "$SAMPLES" | wc -w
qsub $MASTERS_DIR/coverage/slaves/composite.sh
$SAMPLES_DIR $CALL_DIR/samples
invoke file/create.q $DIR $SAMPLES_DIR
qsub $MASTERS_DIR/coverage/slaves/composite_chrom.sh