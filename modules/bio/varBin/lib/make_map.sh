#!/bin/bash

#q    require $SAMPLE $CHROMS $CHROM_SIZES $MIN_MAPQ
#q    require $LIB_DIR $MAPS_DIR $BAM_FILES

#$    -N  vbMap_$SAMPLE
#$    -wd $MAPS_DIR
#$    -l  vf=2G
#$    -t  1-$N_CHROMS

getTaskObject CHROM $CHROMS
getTaskObject CHROM_SIZE $CHROM_SIZES

echo "making map of coverage breaks for $SAMPLE $CHROM"

SAMPLE_DIR=$MAPS_DIR/$SAMPLE
mkdir -p $SAMPLE_DIR
MAP_FILE=$SAMPLE_DIR/$SAMPLE.$CHROM.map
STAT_FILE=$MAP_FILE.stats

slurp samtools merge -u -R $CHROM - $BAM_FILES |
samtools view -q $MIN_MAPQ -f 0x002 -F 0x910 | 
cut -f 4,9 |
awk '$2>0' |
groupBy -g 1 -c 2 -o distinct |
perl $LIB_DIR/make_map.pl $CHROM_SIZE $STAT_FILE |
slurp -o $MAP_FILE
checkPipe

echo
cat $STAT_FILE

echo
echo "done"

#Col	Field	Description
#1	QNAME	Query template/pair NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
#8	MPOS	1-based Mate POSistion
#9	TLEN	inferred Template LENgth (insert size)
#10	SEQ	query SEQuence on the same strand as the reference
#11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
#12+	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE

#Flag	Chr	Description
#0x001	p	the read is paired in sequencing
#0x002	P	the read is mapped in a proper pair
#0x004	u	the query sequence itself is unmapped
#0x008	U	the mate is unmapped
#0x010	r	strand of the query (1 for reverse)
#0x020	R	strand of the mate
#0x040	1	the read is the first read in a pair
#0x080	2	the read is the second read in a pair
#0x100	s	the alignment is not primary
#0x200	f	the read fails platform/vendor quality checks
#0x400	d	the read is either a PCR or an optical duplicate
#0x800      supplementary alignment
