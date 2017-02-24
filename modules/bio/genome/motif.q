
$BIN_SIZE $BIN_SIZE ? $BIN_SIZE : 1000
$NAME     $NAME     ? $NAME     : $MOTIF

invoke bio/genome/file_schema.q
invoke file/require.q $FILE $GENOME_FASTA
$MOTIF_BED  $GENOME_DIR/$GENOME.$NAME.bin_$BIN_SIZE.bed

invoke bio/genome/stat.q $STAT_NAME N_CHROMOSOMES CHROMOSOMES  
qsub lib/motif_array.sh

qsub lib/motif_merge.sh

