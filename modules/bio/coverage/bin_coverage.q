
# collect file structure information
invoke $INSTRS_DIR/file_info/bin_info.q

# calculate bin coverage over all chromosomes
qsub bin_coverage.sh
<file name="bin_coverage.sh">
#!/bin/bash
#q    require $BIN_SAMPLE $BIN_SIZE $NCHROMS $BIN_COV_BAM $BC_BAM $CHROM
#$    -N  binCov_$BIN_SAMPLE\_$BIN_SIZE
#$    -wd $RUN_DIR
#$    -l  vf=5G
#$    -l  h_rt=3:00:00
#$    -t  1-$NCHROMS
#PBS  -N  binCov_$BIN_SAMPLE\_$BIN_SIZE
#PBS  -d  $RUN_DIR
#PBS  -l  mem=5gb
#PBS  -l  walltime=3:00:00
#PBS  -t  1-$NCHROMS
echo "calculating $BIN_SAMPLE $CHROM bin=$BIN_SIZE coverage"
OUT_BAM="$BIN_COV_BAM.chr_"$TASK_ID
checkForData "samtools view $BC_BAM $CHROM"
slurp samtools view -h $BC_BAM $CHROM | 
awk 'BEGIN {
  FS="\t";
  OFS="\t";
  bs='$BIN_SIZE';
}{
  if ($0~/^@/) {             # pass header lines unperturbed
    print $0;
  } else {
    $4=int(($4/bs)+0.5)*bs;  # bin the position
    split($13,fc,":");       # retrieve fractional count for the position
    b=$2"\t"$3"\t"$4;
    fcs[b]+=fc[3];           # sum fractional counts by bin
  }
} END {
  for(b in fcs){             # print and sort the BAM
    print("*",b,255,"1M","*",0,0,"*","*","RG:Z:'$BIN_SAMPLE'","XC:f:"fcs[b]);
  }
}' |
samtools view -hSu - |
samtools sort -o -m 4000000000 - $OUT_BAM.sort |
slurp -o $OUT_BAM
checkPipe
echo "done"
</file>

# merge and index the split bin coverage bams
qsub mrg_idx_bin_cov.sh
<file name="mrg_idx_bin_cov.sh">
#!/bin/bash
#q    require $BIN_SAMPLE $BIN_SIZE $BIN_COV_BAM
#$    -N  mrg_binCov_$BIN_SAMPLE\_$BIN_SIZE
#$    -wd $RUN_DIR
#$    -l  vf=1G
#$    -l  h_rt=2:00:00
#PBS  -N  mrg_binCov_$BIN_SAMPLE\_$BIN_SIZE
#PBS  -d  $RUN_DIR
#PBS  -l  mem=1gb
#PBS  -l  walltime=2:00:00
echo "merging $BIN_SAMPLE bin coverage bam"
BAM_GLOB="$BIN_COV_BAM.chr_*"
samtools merge -f $BIN_COV_BAM $BAM_GLOB
rm $BAM_GLOB
echo "indexing bin coverage bam"
samtools index $BIN_COV_BAM
snipBam $BIN_COV_BAM
echo "done"
</file>


