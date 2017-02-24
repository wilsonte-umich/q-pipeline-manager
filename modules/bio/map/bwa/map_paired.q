
invoke bio/map/bwa/lib/common.q
invoke file/require.q $FILE $GENOME_FASTA $BWT_FILE 

invoke file/require.q $FILE $FASTQ_GLOB1
$N_FILES1 $N_FILES
invoke file/require.q $FILE $FASTQ_GLOB2
$N_FILES2 $N_FILES
$ERROR run [ "$N_FILES1" = "$N_FILES2" ] || \
    echo '"map/bwa/map_paired.q: $FASTQ_GLOB1 and $FASTQ_GLOB2 have different numbers of files"'
dieIf $ERROR
invoke file/create.q $DIR $OUTPUT_DIR

$BAM_FILE $OUTPUT_DIR/$NAME.bam
$BAM_GLOB $BAM_FILE.bwa_*
rm -f $BAM_GLOB 2>/dev/null

$RM_DUP $RM_DUP ? samtools rmdup - - : cat

qsub map.sh

qsub lib/merge.sh

<file name=map.sh>
#!/bin/bash

#q    require $GENOME $GENOME_FASTA $NAME $FASTQ_GLOB1 $FASTQ_GLOB2 $BAM_FILE
#q    option  $ALN_OPTIONS[] $SAMXE_OPTIONS[] $TMP_DIR[/tmp] 

#$    -N  bwa_mapp_$NAME\_$GENOME
#$    -wd $OUTPUT_DIR
#$    -l  vf=4G
#$    -l  h_rt=12:00:00
#$    -t  1-$N_FILES1

#PBS  -N  bwa_mapp_$NAME\_$GENOME
#PBS  -d  $OUTPUT_DIR
#PBS  -l  mem=4gb
#PBS  -l  walltime=12:00:00
#PBS  -t  1-$N_FILES1

getTaskFile $FASTQ_GLOB1
FASTQ_FILE1=$TASK_FILE
getTaskFile $FASTQ_GLOB2
FASTQ_FILE2=$TASK_FILE
BAM_FILE="$BAM_FILE.bwa_$TASK_ID"

echo "running bwa on $NAME, $GENOME, $FASTQ_FILE1 and $FASTQ_FILE2"
mkdir -p $TMP_DIR
TMP_FILE="`echo $BAM_FILE | sed 's|/|_|g'`"

bwa sampe $SAMXE_OPTIONS $GENOME_FASTA \
<(slurp $FASTQ_FILE1 | bwa aln $ALN_OPTIONS $GENOME_FASTA -) \
<(slurp $FASTQ_FILE2 | bwa aln $ALN_OPTIONS $GENOME_FASTA -) \
$FASTQ_FILE1 $FASTQ_FILE2 | 
samtools view -hSb - |
samtools sort -o - $TMP_DIR/$TMP_FILE | 
slurp -o $BAM_FILE

echo "done"

</file>

