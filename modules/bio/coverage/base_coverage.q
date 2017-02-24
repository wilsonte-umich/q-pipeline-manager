
# calculate base coverage per strand
$STRANDED_FLAGS $INVERT_STRANDS ? 0 16 : 16 0  # invert stranded strands since were 1st-strand cDNA
$STRANDS        $STRANDED ? - \+ : $NULL       # escape '+' so that it doesn't get acted on during invoke
$STRAND_FILTERS $STRANDED ? -f16 -F16 : $NULL       
$STRAND_FLAGS   $STRANDED ? $STRANDED_FLAGS : 0        
 
invoke $INSTRS_DIR/count/base_coverage_strand.q \
    $STRAND $STRANDS + $STRAND_FLAG $STRAND_FLAGS + $STRAND_FILTER $STRAND_FILTERS

# merge and index the split base coverage bams
qsub mrg_idx_base_cov.sh
<file name="mrg_idx_base_cov.sh">
#!/bin/bash
#q    require $COUNT_SAMPLE $BC_BAM
#$    -N  mrg_baseCov_$COUNT_SAMPLE
#$    -wd $RUN_DIR
#$    -l  vf=1G
#$    -l  h_rt=3:00:00
#PBS  -N  mrg_baseCov_$COUNT_SAMPLE
#PBS  -d  $RUN_DIR
#PBS  -l  mem=1gb
#PBS  -l  walltime=3:00:00
echo "merging $COUNT_SAMPLE base coverage bam"
BAM_GLOB="$BC_BAM.chr_*_strand_*"
samtools merge -f $BC_BAM $BAM_GLOB
rm $BAM_GLOB
echo "indexing base coverage bam"
samtools index $BC_BAM
snipBam $BC_BAM
echo "done"
</file>


