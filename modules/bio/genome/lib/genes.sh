#!/bin/bash

#q    require $GENOME_DIR $GENOME $ANNOT $TRANS_TXT $GENES_BED $REFGENE_TXT
#q    require $NAMES_TXT $TRANS_ID_COL $GENE_ID_COL

#$    -N  genes_$ANNOT\_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=1G
#$    -l  h_rt=1:00:00

#PBS  -N  genes_$ANNOT\_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=1gb
#PBS  -l  walltime=1:00:00

if [ "$ANNOT" = "UCSC" ]; then
AWK='
g=gs[$1]?gs[$1]:$1;
g=g"::"$2"::"$3;
ss[g]=ss[g]?ss[g]:$4;
ss[g]=ss[g]<=$4?ss[g]:$4;
es[g]=es[g]?es[g]:$5;
es[g]=es[g]>=$5?es[g]:$5;'
else 
AWK='
g=gs[$2]?gs[$2]:$13;
g=g"::"$3"::"$4;
ss[g]=ss[g]?ss[g]:$5;
ss[g]=ss[g]<=$5?ss[g]:$5;
es[g]=es[g]?es[g]:$6;
es[g]=es[g]>=$6?es[g]:$6;'
fi

if [ -f $TRANS_TXT ]; then
    echo "creating $GENOME $ANNOT genes bed file"
    cat <(awk 'BEGIN{FS="\t";OFS="\t"}{print "geneSymbol",$'$TRANS_ID_COL',$'$GENE_ID_COL'}' $NAMES_TXT) $TRANS_TXT |
    awk 'BEGIN{
      FS="\t";
      OFS="\t";
    }{ 
      if($1=="geneSymbol"){ gs[$2]=$3 }else{ '"$AWK"' }
    }END{
      for (g in ss){
        split(g,ga,"::");
        print ga[2], ss[g], es[g], ga[1], 0, ga[3];
      }
    }' |
    sort -k4,4 > $GENES_BED
    checkPipe
    snipFile $GENES_BED
    
    if [ "$TRANS_TXT" = "$REFGENE_TXT" ]; then
        echo "setting genome stats"
        set_stat $GENOME_DIR $GENOME \
          N_GENES,"`cat $GENES_BED | wc -l`" \
          N_GENE_BASES,"`awk '{n+=$3-$2}END{print n?n:0}' $GENES_BED`"
        checkPipe
        snipFile $GENOME_DIR/$GENOME.stats.csv 100 
        echo "done"
    fi
else
    echo "does not exist: $TRANS_TXT"
fi

#-----------------
#UCSC table schema
#-----------------
#REFGENE
#bin	2166	smallint(5) unsigned	range	Indexing field to speed chromosome range queries.
#name	NM_001017364	varchar(255)	values	Name of gene (usually transcript_id from GTF)
#chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
#strand	+	char(1)	values	+ or - for strand
#txStart	207262583	int(10) unsigned	range	Transcription start position
#txEnd	207273337	int(10) unsigned	range	Transcription end position
#cdsStart	207262876	int(10) unsigned	range	Coding region start
#cdsEnd	207273274	int(10) unsigned	range	Coding region end
#exonCount	6	int(10) unsigned	range	Number of exons
#exonStarts	207262583,207263655,2072649...	longblob	 	Exon start positions
#exonEnds	207262934,207263826,2072651...	longblob	 	Exon end positions
#score	0	int(11)	range	score
#name2	C4BPB	varchar(255)	values	Alternate name (e.g. gene_id from GTF)
#cdsStartStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
#cdsEndStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
#exonFrames	0,1,1,1,2,0,	longblob	 	Exon frame {0,1,2}, or -1 if no frame for exon

#ENSEMBL
#bin	585	smallint(5) unsigned	range	Indexing field to speed chromosome range queries.
#name	ENST00000456328	varchar(255)	values	Name of gene (usually transcript_id from GTF)
#chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
#strand	+	char(1)	values	+ or - for strand
#txStart	11868	int(10) unsigned	range	Transcription start position
#txEnd	14409	int(10) unsigned	range	Transcription end position
#cdsStart	14409	int(10) unsigned	range	Coding region start
#cdsEnd	14409	int(10) unsigned	range	Coding region end
#exonCount	3	int(10) unsigned	range	Number of exons
#exonStarts	11868,12612,13220,	longblob	 	Exon start positions
#exonEnds	12227,12721,14409,	longblob	 	Exon end positions
#score	0	int(11)	range	score
#name2	ENSG00000223972	varchar(255)	values	Alternate name (e.g. gene_id from GTF)
#cdsStartStat	none	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
#cdsEndStat	none	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
#exonFrames	-1,-1,-1,	longblob	 	Exon frame {0,1,2}, or -1 if no frame for exon

#UCSC
#name	uc001aaa.3	varchar(255)	values	Name of gene
#chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
#strand	+	char(1)	values	+ or - for strand
#txStart	11873	int(10) unsigned	range	Transcription start position
#txEnd	14409	int(10) unsigned	range	Transcription end position
#cdsStart	11873	int(10) unsigned	range	Coding region start
#cdsEnd	11873	int(10) unsigned	range	Coding region end
#exonCount	3	int(10) unsigned	range	Number of exons
#exonStarts	11873,12612,13220,	longblob	 	Exon start positions
#exonEnds	12227,12721,14409,	longblob	 	Exon end positions
#proteinID	 	varchar(40)	values	UniProt display ID for Known Genes, UniProt accession or RefSeq protein ID for UCSC Genes
#alignID	uc001aaa.3	varchar(255)	values	Unique identifier for each (known gene, alignment position) pair

