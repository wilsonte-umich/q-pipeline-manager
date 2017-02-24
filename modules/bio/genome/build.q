
echo "------------------------------"
echo $NAME
echo "------------------------------"

# input genome information
invoke bio/genome/file_schema.q
invoke bio/genome/stat.q $STAT_NAME CHROMOSOMES N_CHROMOSOMES SEQUENCES N_SEQUENCES 
$INPUT_GENOME     $GENOME
$INPUT_GENOME_DIR $GENOME_DIR
invoke file/require.q $FILE $INPUT_GENOME_DIR/chromosomes/*.fa

# bed files to renumber and populate
>>>
$ALWAYS_LIFT 
run ls -1 
$INPUT_GENOME_DIR/$INPUT_GENOME.gap.bed
$INPUT_GENOME_DIR/$INPUT_GENOME.*.transcripts.bed
$INPUT_GENOME_DIR/$INPUT_GENOME.*.genes.bed
$INPUT_GENOME_DIR/$INPUT_GENOME.*.map.*.bed
2>/dev/null    
<<<
$LIFT_BEDS $LIFT_BEDS ? $LIFT_BEDS $ALWAYS_LIFT : $ALWAYS_LIFT

# output genome information
$OUTPUT_GENOME     $NAME
$OUTPUT_GENOME_DIR $GENOMES_DIR/$OUTPUT_GENOME
invoke file/create.q $DIR $OUTPUT_GENOME_DIR/chromosomes
$TARGET_DIR   $OUTPUT_GENOME_DIR
$FILE_ROOT    $OUTPUT_GENOME
$STAT_NAME    GENOME
$STAT_VALUE   $OUTPUT_GENOME
invoke stat/set.q  
$STAT_NAME    PARENT_GENOME
$STAT_VALUE   $INPUT_GENOME
invoke stat/set.q 
$STAT_NAME    CHROMOSOMES
$STAT_VALUE   $CHROMOSOMES
invoke stat/set.q  
$STAT_NAME    N_CHROMOSOMES
$STAT_VALUE   $N_CHROMOSOMES
invoke stat/set.q  
$STAT_NAME    SEQUENCES
$STAT_VALUE   $SEQUENCES
invoke stat/set.q  
$STAT_NAME    N_SEQUENCES
$STAT_VALUE   $N_SEQUENCES
invoke stat/set.q  

# change information
invoke file/require.q $FILE $CHANGES_FILE
$FEATURES_FILE $FEATURES_FILE ? $FEATURES_FILE : $NULL

# modify chromosomes
qsub lib/build_chrom.pl $CHROMOSOME $CHROMOSOMES

# merge chromosomes
qsub lib/build_merge.pl

