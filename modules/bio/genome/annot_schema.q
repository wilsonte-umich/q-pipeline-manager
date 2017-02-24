
# files containing the RefSeq gene annotation
$TRANSCRIPTS_TXT  $GENOME_DIR/$ANNOTATION.txt  # as obtained from UCSC; transcripts
$TRANSCRIPTS_BED  $GENOME_DIR/$GENOME.$ANNOTATION.transcripts.bed   # transcript spans from $TRANSCRIPTS_TXT, one line per transcript
$GENES_BED        $GENOME_DIR/$GENOME.$ANNOTATION.genes.bed         # gene spans from $TRANSCRIPTS_TXT, one line per gene (i.e. all of its transcripts)
$MAP_STRANDED     $GENOME_DIR/$GENOME.$ANNOTATION.map.stranded.bed  # exons, introns, etc.
$MAP_UNSTRANDED   $GENOME_DIR/$GENOME.$ANNOTATION.map.unstranded.bed  
$MAP_GENES_BED    $GENOME_DIR/$GENOME.$ANNOTATION.map.genes.bed     # gene spans from $MAP_STRANDED

preserve all

