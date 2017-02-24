#!/bin/bash

#q    require $GENOME_DIR $GENOME $GAP_TXT $GAP_BED

#$    -N  gaps_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=1G
#$    -l  h_rt=1:00:00

#PBS  -N  gaps_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=1gb
#PBS  -l  walltime=1:00:00

if [ -f $GAP_TXT ]; then

    echo "creating $GENOME gap bed file"
    awk 'BEGIN{
      FS="\t";
      OFS="\t";
    }{
      print $2, $3, $4, $8, $7, "+";
    }' $GAP_TXT > $GAP_BED
    checkPipe
    snipFile $GAP_BED
    
    echo "setting genome stats"
    set_stat $GENOME_DIR $GENOME \
      N_GAPS,"`cat $GAP_BED | wc -l`" \
      N_GAP_BASES,"`awk '{n+=$3-$2}END{print n?n:0}' $GAP_BED`"
    snipFile $GENOME_DIR/$GENOME.stats.csv 100
    checkPipe
    echo "done"
    
else
    echo "does not exist: $GAP_TXT"
fi

#-----------------
#UCSC table schema
#-----------------
#field	example	SQL type	info	description
#bin	585	smallint(6)	range	Indexing field to speed chromosome range queries.
#chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
#chromStart	0	int(10) unsigned	range	start position in chromosome
#chromEnd	10000	int(10) unsigned	range	end position in chromosome
#ix	1	int(11)	range	ix of this fragment (useless)
#n	N	char(1)	values	always 'N'
#size	10000	int(10) unsigned	range	size of gap
#type	telomere	varchar(255)	values	contig, clone, fragment, etc.
#bridge	no	varchar(255)	values	yes, no, mrna, bacEndPair, etc.

