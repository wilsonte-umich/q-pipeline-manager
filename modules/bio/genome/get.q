
echo "------------------------------"
echo $GENOME
echo "------------------------------"

invoke bio/genome/file_schema.q
invoke file/create.q $DIR $GENOME_DIR

$GOLDEN_PATH  ftp://hgdownload.cse.ucsc.edu/goldenPath

thread get
$TYPE     bigZip
$URL      $GOLDEN_PATH/$GENOME/bigZips
$BZ_FILES $BZ_FILES ? $BZ_FILES : $NULL
$FILES \
    chromFa.tar.gz \
    chromFa.zip \
    refMrna.fa.gz \
    $BZ_FILES
qsub lib/get.sh
  
# do not thread: UCSC will balk at too many sessions  
$TYPE     annotation  
$URL      $GOLDEN_PATH/$GENOME/database
$DB_FILES $DB_FILES ? $DB_FILES : $NULL
$FILES \
    chromInfo.txt.gz \
    gap.txt.gz \
    refGene.txt.gz \
    knownGene.txt.gz \
    ensGene.txt.gz \
    ensemblToGeneName.txt.gz \
    kgXref.txt.gz \
    $DB_FILES
qsub lib/get.sh

$ANNOTS          refGene                  UCSC                   Ensembl
$NAMES_TXTS      $UCSC_NAMES_TXT          $UCSC_NAMES_TXT        $ENSEMBL_NAMES_TXT
$TRANS_ID_COLS   6                        1                      1
$GENE_ID_COLS    5                        5                      2
$TRANS_TXTS      $REFGENE_TXT             $UCSC_TXT              $ENSEMBL_TXT
$GENES_BEDS      $REFGENE_GENES_BED       $UCSC_GENES_BED        $ENSEMBL_GENES_BED
$TRANS_BEDS      $REFGENE_TRANSCRIPTS_BED $UCSC_TRANSCRIPTS_BED  $ENSEMBL_TRANSCRIPTS_BED
$MAP_STR_BEDS    $REFGENE_MAP_STRANDED    $UCSC_MAP_STRANDED     $ENSEMBL_MAP_STRANDED
$MAP_UNSTR_BEDS  $REFGENE_MAP_UNSTRANDED  $UCSC_MAP_UNSTRANDED   $ENSEMBL_MAP_UNSTRANDED
$MAP_GENES_BEDS  $REFGENE_MAP_GENES_BED   $UCSC_MAP_GENES_BED    $ENSEMBL_MAP_GENES_BED

thread parse get
qsub lib/fasta.sh
qsub lib/gaps.sh
qsub lib/genes.sh       $ANNOT $ANNOTS + $TRANS_TXT $TRANS_TXTS + $GENES_BED $GENES_BEDS + $NAMES_TXT $NAMES_TXTS + $TRANS_ID_COL $TRANS_ID_COLS + $GENE_ID_COL $GENE_ID_COLS
qsub lib/transcripts.sh $ANNOT $ANNOTS + $TRANS_TXT $TRANS_TXTS + $TRANS_BED $TRANS_BEDS + $NAMES_TXT $NAMES_TXTS + $TRANS_ID_COL $TRANS_ID_COLS + $GENE_ID_COL $GENE_ID_COLS
qsub lib/chunks.sh

thread map get
qsub lib/map.pl  $ANNOT $ANNOTS + $TRANS_TXT $TRANS_TXTS + $MAP_STR_BED $MAP_STR_BEDS + $MAP_UNSTR_BED $MAP_UNSTR_BEDS + $MAP_GENES_BED $MAP_GENES_BEDS + $NAMES_TXT $NAMES_TXTS + $TRANS_ID_COL $TRANS_ID_COLS + $GENE_ID_COL $GENE_ID_COLS

