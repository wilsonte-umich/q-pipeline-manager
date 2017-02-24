#!/bin/bash

#q    require $GENOME $NAME $BIN_SIZE $MOTIF $CHROMOSOMES $GENOME_FASTA $MOTIF_BED

#$    -N  $NAME\_$BIN_SIZE\_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=2G
#$    -l  h_rt=2:00:00
#$    -t  1-$N_CHROMOSOMES

#PBS  -N  $NAME\_$BIN_SIZE\_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=2gb
#PBS  -l  walltime=2:00:00
#PBS  -t  1-$N_CHROMOSOMES

getTaskObject CHROMOSOME $CHROMOSOMES
echo "creating BED file of $GENOME $CHROMOSOME $NAME content at bin size $BIN_SIZE"    
OUT_FILE=$MOTIF_BED.chrom_$CHROMOSOME

slurp samtools faidx $GENOME_FASTA $CHROMOSOME | 
grep -vP '^>' | 
awk '{printf $0}' |
perl -e '
    $m = "'$MOTIF'";    
    $bs = '$BIN_SIZE';
    $wbin = 0;
    $rf = 10000;
    $tmp = $m;
    while($tmp =~ s~\w\|\w~x~g) {}
    $tmp =~ s~[\(|\)]~~g;
    $nPre = length($tmp) - 1;
    $pre = "";        
    while(read STDIN, $seq, $bs){
        $seq = uc($seq);
        $tgt = "$pre$seq";
        $val = ($tgt =~ s/$m//g);
        $val or $val = 0;
        $val = $val / $bs;
        $val = int($val*$rf+0.5)/$rf;
        print join("\t", "'$CHROMOSOME'", $wbin, $wbin + $bs, "*", $val, "+"), "\n";
        $wbin += $bs;
        $nPre and $pre = substr($seq, -$nPre);
    }
' | 
slurp -o $OUT_FILE
checkPipe
snipFile $OUT_FILE

echo "done"

