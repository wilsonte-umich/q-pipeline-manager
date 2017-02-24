
# the directory that contains the files for a specific genome
$GENOME_DIR      $GENOMES_DIR/$GENOME 

# FASTA file containing all and only canonical chromosome sequences
$GENOME_FASTA    $GENOME_DIR/$GENOME.fa

# tab-delimited table with fields: chromosome, size, 2-bit-file
$CHROM_INFO_TXT  $GENOME_DIR/chromInfo.txt  # as obtained from UCSC
$CHROM_SIZES_TXT $GENOME_DIR/$GENOME.chromSizes.txt  # only canonical chromosomes

# files containing genome gaps, i.e. stretches of N bases
$GAP_TXT         $GENOME_DIR/gap.txt  # as obtained from UCSC
$GAP_BED         $GENOME_DIR/$GENOME.gap.bed  # BED format

# files containing the RefSeq gene annotation
$REF_MRNA_FASTA          $GENOME_DIR/refMrna.fa   # as obtained from UCSC
$REFGENE_TXT             $GENOME_DIR/refGene.txt  # as obtained from UCSC; transcripts
$REFGENE_TRANSCRIPTS_BED $GENOME_DIR/$GENOME.refGene.transcripts.bed   # transcript spans from $REFGENE_TXT, one line per transcript
$REFGENE_GENES_BED       $GENOME_DIR/$GENOME.refGene.genes.bed         # gene spans from $REFGENE_TXT, one line per gene (i.e. all of its transcripts)
$REFGENE_MAP_STRANDED    $GENOME_DIR/$GENOME.refGene.map.stranded.bed  # exons, introns, etc.
$REFGENE_MAP_UNSTRANDED  $GENOME_DIR/$GENOME.refGene.map.unstranded.bed  
$REFGENE_MAP_GENES_BED   $GENOME_DIR/$GENOME.refGene.map.genes.bed     # gene spans from $REFGENE_MAP_STRANDED

# files containing the UCSC gene annotation
$UCSC_TXT             $GENOME_DIR/knownGene.txt  # as obtained from UCSC; transcripts
$UCSC_NAMES_TXT       $GENOME_DIR/kgXref.txt     # name conversion file (1=kgID=UCSC ID, 6=refseq=RefSeqID, 5=geneSymbol)
$UCSC_TRANSCRIPTS_BED $GENOME_DIR/$GENOME.UCSC.transcripts.bed   # transcript spans from $UCSC_TXT, one line per transcript
$UCSC_GENES_BED       $GENOME_DIR/$GENOME.UCSC.genes.bed         # gene spans from $UCSC_TXT, one line per gene (i.e. all of its transcripts)
$UCSC_MAP_STRANDED    $GENOME_DIR/$GENOME.UCSC.map.stranded.bed  # exons, introns, etc.
$UCSC_MAP_UNSTRANDED  $GENOME_DIR/$GENOME.UCSC.map.unstranded.bed  
$UCSC_MAP_GENES_BED   $GENOME_DIR/$GENOME.UCSC.map.genes.bed     # gene spans from $UCSC_MAP_STRANDED

# files containing the Ensembl gene annotation
$ENSEMBL_TXT             $GENOME_DIR/ensGene.txt            # as obtained from UCSC; transcripts
$ENSEMBL_NAMES_TXT       $GENOME_DIR/ensemblToGeneName.txt  # name conversion file (1=name=ENST transcriptID, 2=value=geneSymbol)
$ENSEMBL_TRANSCRIPTS_BED $GENOME_DIR/$GENOME.ensembl.transcripts.bed   # transcript spans from $ENSEMBL_TXT, one line per transcript
$ENSEMBL_GENES_BED       $GENOME_DIR/$GENOME.ensembl.genes.bed         # gene spans from $ENSEMBL_TXT, one line per gene (i.e. all of its transcripts)
$ENSEMBL_MAP_STRANDED    $GENOME_DIR/$GENOME.ensembl.map.stranded.bed  # exons, introns, etc.
$ENSEMBL_MAP_UNSTRANDED  $GENOME_DIR/$GENOME.ensembl.map.unstranded.bed  
$ENSEMBL_MAP_GENES_BED   $GENOME_DIR/$GENOME.ensembl.map.genes.bed     # gene spans from $ENSEMBL_MAP_STRANDED

# files containing information on genome chunks
$CHUNKS_BED      $GENOME_DIR/$GENOME.chunks.bed  # BED format

preserve all

