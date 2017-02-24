
invoke bio/map/bwa/lib/common.q
invoke file/require.q $FILE $GENOME_FASTA $BWT_FILE 

invoke file/require.q $FILE $FASTQ_GLOB1
$N_FILES1 $N_FILES
invoke file/create.q $DIR $OUTPUT_DIR

$BAM_FILE $OUTPUT_DIR/$NAME.bam
$BAM_GLOB $BAM_FILE.bwa_*
rm -f $BAM_GLOB 2>/dev/null

$RM_DUP $RM_DUP ? samtools rmdup -s - - : cat

qsub map.sh

qsub lib/merge.sh

<file name=map.sh>
#!/bin/bash

#q    require $GENOME $GENOME_FASTA $NAME $FASTQ_GLOB1 $BAM_FILE
#q    option  $ALN_OPTIONS[] $SAMXE_OPTIONS[] $TMP_DIR[/tmp] 

#$    -N  bwa_maps_$NAME\_$GENOME
#$    -wd $OUTPUT_DIR
#$    -l  vf=4G
#$    -l  h_rt=12:00:00
#$    -t  1-$N_FILES1

#PBS  -N  bwa_maps_$NAME\_$GENOME
#PBS  -d  $OUTPUT_DIR
#PBS  -l  mem=4gb
#PBS  -l  walltime=12:00:00
#PBS  -t  1-$N_FILES1

getTaskFile $FASTQ_GLOB1
FASTQ_FILE1=$TASK_FILE
BAM_FILE="$BAM_FILE.bwa_$TASK_ID"

echo "running bwa on $NAME, $GENOME, $FASTQ_FILE1"
mkdir -p $TMP_DIR
TMP_FILE="`echo $BAM_FILE | sed 's|/|_|g'`"

slurp $FASTQ_FILE1 | 
bwa aln $ALN_OPTIONS $GENOME_FASTA - | 
bwa samse $SAMXE_OPTIONS $GENOME_FASTA - $FASTQ_FILE1 |
samtools view -hSb - |
samtools sort -o - $TMP_DIR/$TMP_FILE | 
slurp -o $BAM_FILE

echo "done"

</file>

