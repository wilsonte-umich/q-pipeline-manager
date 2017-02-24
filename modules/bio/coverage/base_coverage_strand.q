
$STRAND_OPTION     $STRAND ? -strand $STRAND : $NULL

qsub base_coverage_strand.sh

<file name="base_coverage_strand.sh">
#!/bin/bash
#q    require $COUNT_SAMPLE $STRAND $NCHROMS $BC_BAM $STRAND_FILTER $MAP_BAM $CHROM $STRAND_OPTION $EFF_LENGTH $STRAND_FLAG
#$    -N  baseCov_$COUNT_SAMPLE\_$STRAND
#$    -wd $RUN_DIR
#$    -l  vf=4G
#$    -l  h_rt=3:00:00
#$    -t  1-$NCHROMS
#PBS  -N  baseCov_$COUNT_SAMPLE\_$STRAND
#PBS  -d  $RUN_DIR
#PBS  -l  mem=4gb
#PBS  -l  walltime=3:00:00
#PBS  -t  1-$NCHROMS

checkTaskID

echo "creating $NAME $CHROM$STRAND bin coverage bam"

OUT_BAM="$BC_BAM.chr_"$TASK_ID"_strand_$STRAND"

checkForData "samtools view $STRAND_FILTER $MAP_BAM $CHROM"

slurp samtools view -hu $STRAND_FILTER $BAM_FILE $CHROM |

genomeCoverageBed -bg -split $STRAND_OPTION -ibam stdin | # fastest genomeCoverageBed, bedGraph format

awk -v header="`samtools view -H $BAM_FILE`" '            # convert bedGraph to SAM with line for every base position hit
BEGIN {
  FS="\t";
  OFS="\t";
  print header;
}{
  fc=$4/'$READ_LENGTH';                                   # fractionally spread single read count over all its bases   
  for(pos=$2+1;pos<=$3;pos++){                            # include fractional count as tag XC:f: in base coverage bam
    print("*",'$STRAND_FLAG',$1,pos,255,"1M","*",0,0,"*","*","RG:Z:'$COUNT_SAMPLE'","XC:f:"fc);
  }
}' |
samtools view -hSb - |  # parsed bedGraph output is already sorted per chromosome file
slurp -o $OUT_BAM
checkPipe
echo "done"
</file>


