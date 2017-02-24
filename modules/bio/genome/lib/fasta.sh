#!/bin/bash

#q    require $GENOME_DIR $GENOME $GENOME_FASTA $REF_MRNA_FASTA

#$    -N  fasta_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=2G
#$    -l  h_rt=1:00:00

#PBS  -N  fasta_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=2gb
#PBS  -l  walltime=1:00:00

N_CHROM_FILES="`ls -1 $GENOME_DIR/chr*.fa | wc -l`"
if [ "$N_CHROM_FILES" = "0" ]; then
    echo "did not find any chr*.fa files in $GENOME_DIR"
    echo "check the get_$GENOME""_bigZip log file for download errors"
    exit 100
fi

echo "moving individual chromosome files"
CHROM_DIR="$GENOME_DIR/chromosomes"
mkdir -p $CHROM_DIR
mv $GENOME_DIR/chr*.fa $CHROM_DIR
checkPipe

echo "creating and indexing $GENOME_FASTA"
cd $CHROM_DIR
FILES="`ls -1 chr*.fa | grep -v '_'`"
CHROMOSOMES="`ls -1 chr*.fa | grep -v '_' | sed 's|.fa$||'`"
SEQUENCES="`ls -1 *.fa | sed 's|.fa$||'`"
cat $FILES > $GENOME_FASTA
checkPipe
cd $GENOME_DIR
samtools faidx $GENOME_FASTA
checkPipe
echo
echo $GENOME_FASTA
head $GENOME_FASTA
echo "..."
echo

if [ -f $REF_MRNA_FASTA ]; then
    echo "indexing $REF_MRNA_FASTA"
    samtools faidx $REF_MRNA_FASTA
    checkPipe
    N_TRANSCRIPTS="`slurp $REF_MRNA_FASTA | grep -P '^>' | wc -l`"
    echo "$N_TRANSCRIPTS transcripts"
    snipFile $REF_MRNA_FASTA
else
    echo "does not exist: $REF_MRNA_FASTA"
    echo
fi

echo "setting genome stats"
set_stat $GENOME_DIR $GENOME \
  GENOME,$GENOME \
  CHROMOSOMES,"$CHROMOSOMES" \
  N_CHROMOSOMES,"`echo $CHROMOSOMES | wc -w`" \
  SEQUENCES,"$SEQUENCES" \
  N_SEQUENCES,"`echo $SEQUENCES | wc -w`" \
  N_BASES,"`slurp $GENOME_FASTA | grep -vP '^>' | sed 's/\s//g' | awk '{c+=length($0)}END{print c}'`" \
  N_TRANSCRIPTS,$N_TRANSCRIPTS
checkPipe
snipFile $GENOME_DIR/$GENOME.stats.csv 100

echo "done"

