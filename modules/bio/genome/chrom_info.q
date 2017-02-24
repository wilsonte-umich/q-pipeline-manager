
$GENOME_DIR     $GENOMES_DIR/$GENOME
$GENOME_FASTA   $GENOME_DIR/$GENOME.fa
$CHROM_FILTER   awk '\$1!="chrM"&&\$1!~/_/' | sort -k1,1
$CHROM_FILE     $GENOME_DIR/chromInfo.txt
$CHROMS         run cat $CHROM_FILE | $CHROM_FILTER | cut -f1
$CHROM_SIZES    run cat $CHROM_FILE | $CHROM_FILTER | cut -f2
$GAP_FILE       $GENOME_DIR/$GENOME.gap.bed 
$GAP_SIZES      awk -F "\\t" '{cs[\$1]+=\$3-\$2}END{for(c in cs){print c"\\t"cs[c]}}'
$GAP_SIZES      run cat $GAP_FILE | $GAP_SIZES | $CHROM_FILTER | cut -f2 
$GENOME_SIZE    run perl -e 'map{\$gs+=\$_}split(" ","$CHROM_SIZES");print \$gs'
$GAP_SIZE       run perl -e 'map{\$gs+=\$_}split(" ","$GAP_SIZES");print \$gs'
$NOGAP_SIZE     run echo \$(($GENOME_SIZE - $GAP_SIZE))
$N_CHROMS       run echo "$CHROMS" | wc -w
$MAX_CHROM_SIZE run perl -e '(\$m)=sort{\$b<=>\$a}split(" ","$CHROM_SIZES");print \$m'

preserve all
