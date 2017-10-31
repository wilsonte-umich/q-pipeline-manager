
echo "------------------------------"
echo $NAME
echo "------------------------------"

# input genome information
invoke bio/genome/file_schema.q
invoke bio/genome/stat.q $STAT_NAME CHROMOSOMES 
invoke file/require.q $FILE $GENOME_DIR/chromosomes/*.fa
$SPLIT_PREFIX $GENOMES_DIR/$GENOME/$GENOME.split.$NAME

# output information
$IS_FASTA run [ "$OUTPUT_TYPE" = "FASTA" ] && echo 1
$SUFFIX   $IS_FASTA ? fa : bed
$IS_FASTA $IS_FASTA ? $IS_FASTA : 0
$IS_BED   run [ "$OUTPUT_TYPE" = "BED" ] && echo 1
$IS_BED   $IS_BED ? $IS_BED : 0
$CHECK_SUM run echo \$(($IS_FASTA + $IS_BED))
dieIf    run [ "$CHECK_SUM" = "0" ] && echo "Unrecognized \$OUTPUT_TYPE"
$OUT_FILE $SPLIT_PREFIX.$SUFFIX 

# process chromosomes
qsub lib/split_chrom.pl $CHROMOSOME $CHROMOSOMES

# merge chromosomes
qsub lib/split_merge.pl

