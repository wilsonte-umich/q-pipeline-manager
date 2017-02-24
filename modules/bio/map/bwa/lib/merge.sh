#!/bin/bash

#q    require $NAME $GENOME $BAM_FILE $N_FILES1 $BAM_GLOB $RM_DUP

#$    -N  bwa_mrg_$NAME\_$GENOME
#$    -wd $OUTPUT_DIR
#$    -l  vf=1G
#$    -l  h_rt=6:00:00

#PBS  -N  bwa_mrg_$NAME\_$GENOME
#PBS  -d  $OUTPUT_DIR
#PBS  -l  mem=1gb
#PBS  -l  walltime=6:00:00

if [ "$N_FILES1" = "1" ]; then
    IN_BAM="$BAM_FILE.bwa_1"
    if [ "$RM_DUP" = "cat" ]; then
        mv $IN_BAM $BAM_FILE
    else 
        slurp $IN_BAM |
        $RM_DUP |
        slurp -o $BAM_FILE
        rm $IN_BAM
    fi
else 
    echo "merging $NAME bwa bam files"
    samtools merge - $BAM_GLOB |
    $RM_DUP |
    slurp -o $BAM_FILE
    rm $BAM_GLOB
fi

echo "indexing $BAM_FILE"
samtools index $BAM_FILE

echo
echo "$BAM_FILE"
samtools view -H $BAM_FILE
samtools view $BAM_FILE | head
echo ...
echo

echo "done"

