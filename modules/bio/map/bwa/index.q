
invoke bio/map/bwa/lib/common.q
invoke file/require.q $FILE $GENOME_FASTA

qsub index.sh

<file name=index.sh>
#!/bin/bash

#q    require $GENOME $GENOME_FASTA

#$    -N  bwa_index_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=4G
#$    -l  h_rt=12:00:00

#PBS  -N  bwa_index_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=4gb
#PBS  -l  walltime=12:00:00

echo "creating bwa index for genome $GENOME"
snipFile $GENOME_FASTA
bwa index $GENOME_FASTA
echo "done"

</file>

