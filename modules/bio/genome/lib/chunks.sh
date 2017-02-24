#!/bin/bash

#q    require $GENOME_DIR $GENOME $CHROM_INFO_TXT $CHROM_SIZES_TXT $CHUNKS_BED
#q    option  $N_CHUNKS[75]

#$    -N  chunks_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=1G
#$    -l  h_rt=1:00:00

#PBS  -N  chunks_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=1gb
#PBS  -l  walltime=1:00:00

if [ -f $CHROM_INFO_TXT ]; then

    echo "writing file of canonical chromosome sizes"
    cat $CHROM_INFO_TXT | 
    cut -f1,2 | 
    grep -P '^chr' | 
    grep -v '_' > $CHROM_SIZES_TXT  
    checkPipe
    snipFile $CHROM_SIZES_TXT
    
    N_BASES="`get_stat $GENOME_DIR $GENOME N_BASES`"      
    if [ "$N_BASES" = "" ] || [ "$N_BASES" = "0" ]; then
        echo "failed to recover genome length from stats file"
        echo "check the get_$GENOME""_bigZip log file for download errors"
        exit 100
    fi

    echo "splitting $GENOME into ~$N_CHUNKS chunks"
    CHUNK_SIZE="`perl -e '
      $x=int('$N_BASES'/'$N_CHUNKS'); 
      $s=10**(length($x)-2); 
      print int($x/$s+0.5)*$s';
    `"
    cat $CHROM_SIZES_TXT | 
    awk '{ 
      z='$CHUNK_SIZE';
      for(s=1;s<=$2;s+=z){
        r=$2-s+1;
        e=r>z*1.1?s+z-1:$2;
        print $1"\t"s-1"\t"e
        if(e==$2)break;
      }
    }' > $CHUNKS_BED
    checkPipe
    snipFile $CHUNKS_BED

    echo "setting genome stats"
    set_stat $GENOME_DIR $GENOME \
      NOMINAL_CHUNK_SIZE,$CHUNK_SIZE \
      N_CHUNKS,"`cat $CHUNKS_BED | wc -l`" \
      CHUNKS,"`awk '{print $1":"$2+1"-"$3}' $CHUNKS_BED`" \
      CHUNK_CHROMOSOMES,"`awk '{print $1}' $CHUNKS_BED`" \
      CHUNK_STARTS,"`awk '{print $2+1}' $CHUNKS_BED`" \
      CHUNK_ENDS,"`awk '{print $3}' $CHUNKS_BED`" \
      CHUNK_SIZES,"`awk '{print $3-$2}' $CHUNKS_BED`" 
    checkPipe
    snipFile $GENOME_DIR/$GENOME.stats.csv 100
    
else
    echo "does not exist: $CHROM_INFO_TXT"
fi

#-----------------
#UCSC table schema
#-----------------
#field	example	SQL type	info	description
#bin	790	smallint(5) unsigned	range	Indexing field to speed chromosome range queries.
#name	NR_026775	varchar(255)	values	Name of gene (usually transcript_id from GTF)
#chrom	chr6	varchar(255)	values	Reference sequence chromosome or scaffold
#strand	+	char(1)	values	+ or - for strand
#txStart	26924771	int(10) unsigned	range	Transcription start position
#txEnd	26991753	int(10) unsigned	range	Transcription end position
#cdsStart	26991753	int(10) unsigned	range	Coding region start
#cdsEnd	26991753	int(10) unsigned	range	Coding region end
#exonCount	3	int(10) unsigned	range	Number of exons
#exonStarts	26924771,26936236,26991032,	longblob	 	Exon start positions
#exonEnds	26924962,26936301,26991753,	longblob	 	Exon end positions
#score	0	int(11)	range	score
#name2	LINC00240	varchar(255)	values	Alternate name (e.g. gene_id from GTF)
#cdsStartStat	unk	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
#cdsEndStat	unk	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
#exonFrames	-1,-1,-1,	longblob	 	Exon frame {0,1,2}, or -1 if no frame for exon

