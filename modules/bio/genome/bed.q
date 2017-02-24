
invoke bio/genome/file_schema.q
invoke file/create.q $DIR $GENOME_DIR

$GOLDEN_PATH  ftp://hgdownload.cse.ucsc.edu/goldenPath
$TYPE   $NAME  
$URL    $GOLDEN_PATH/$GENOME/database
$FILES  $UCSC_FILE
qsub lib/get.sh

$UCSC_FILE $GENOME_DIR/$UCSC_FILE
$BED_FILE  $GENOME_DIR/$GENOME.$NAME.bed
qsub lib/bed.sh

