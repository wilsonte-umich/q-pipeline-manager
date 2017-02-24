
invoke bio/genome/file_schema.q
invoke bio/genome/stat.q $STAT_NAME N_CHUNKS CHUNKS CHUNK_SIZES
qsub lib/chunk_read_coverage.sh
qsub lib/merge_read_coverage.sh

#invoke lib/merge_read_coverage.q

## merge and index the split base coverage bams
#qsub mrg_idx_base_cov.sh
#<file name="mrg_idx_base_cov.sh">
##!/bin/bash
##q    require $COUNT_SAMPLE $BC_BAM
##$    -N  mrg_baseCov_$COUNT_SAMPLE
##$    -wd $RUN_DIR
##$    -l  vf=1G
##$    -l  h_rt=3:00:00
##PBS  -N  mrg_baseCov_$COUNT_SAMPLE
##PBS  -d  $RUN_DIR
##PBS  -l  mem=1gb
##PBS  -l  walltime=3:00:00
#echo "merging $COUNT_SAMPLE base coverage bam"
#BAM_GLOB="$BC_BAM.chr_*_strand_*"
#samtools merge -f $BC_BAM $BAM_GLOB
#rm $BAM_GLOB
#echo "indexing base coverage bam"
#samtools index $BC_BAM
#snipBam $BC_BAM
#echo "done"
#</file>


