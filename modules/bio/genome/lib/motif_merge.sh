#!/bin/bash

#q    require $GENOME $NAME $BIN_SIZE $MOTIF_BED

#$    -N  mrg_$NAME\_$BIN_SIZE\_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=2G
#$    -l  h_rt=2:00:00

#PBS  -N  mrg_$NAME\_$BIN_SIZE\_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=2gb
#PBS  -l  walltime=2:00:00

echo "merging BED file of $GENOME $NAME content at bin size $BIN_SIZE"   

IN_GLOB=$MOTIF_BED.chrom_*
slurp -o $MOTIF_BED $IN_GLOB
rm $IN_GLOB
snipFile $MOTIF_BED

echo "done"

