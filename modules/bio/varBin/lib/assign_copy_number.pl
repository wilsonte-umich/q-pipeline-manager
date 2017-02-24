use strict;
use warnings;

# using bin widths and counts to determine copy numbers by Hidden Markov

# require common subs
require "$ENV{LIB_DIR}/common.pl";
use vars qw(@bins %binFileCols $nFileLeadCol);

# global variables
our $GROUP_NAME = $ENV{GROUP_NAME};
our @SAMPLES = split(" ", $ENV{SAMPLES});
our @CHROMS  = split(" ", $ENV{CHROMS});
our $MAX_CN  = $ENV{MAX_CN};

# get the reference chromosome information
print STDERR "loading reference statistics\n";
getStats($ENV{REF_STAT_FILE}, \our %refStats);
our $REF_PLOIDY  = $refStats{REF_PLOIDY}[0];
my $refMean     = $refStats{refMean}[0];
my $refStdev    = $refStats{refStdev}[0];

# set the boundaries for the bin-width HMM
print STDERR "setting HMM boundary limits\n";
my $nSD   = 3.5;  # helps determine HMM boundary limits
my $nBins = 100; # number of divisions, i.e. observation states, in the HMM (actually get one more)
my $maxBinI = $nBins + 2;
my ($mean1, $stdev1) = adjustSizeStats($refMean, $refStdev, 1);       # CN=1, largest varBin size
my ($meanX, $stdevX) = adjustSizeStats($refMean, $refStdev, $MAX_CN); # smallest varBin size
my $maxSize = $mean1 + $nSD * $stdev1; # largest and smallest varBin sizes to allow in HMM
my $minSize = $meanX - $nSD * $stdevX;
my $binWidth = int(($maxSize - $minSize) / $nBins + 0.5); # width of each bin in unit varBin size
my $maxBinSize = getBinSize($maxSize); # rounded bin values at the limit bins in unit varBin length
my $minBinSize = getBinSize($minSize);
print STDERR "    CN:1-$MAX_CN, bins:$minBinSize-$maxBinSize($binWidth)\n";

# generate a set of copy number emission probabilities
print STDERR "generating HMM emission probabilities\n";
my $eProbFile    = "$ENV{BINS_DIR}/$GROUP_NAME.emiss.prob.txt";
$ENV{eProbFile}  = $eProbFile;
$ENV{refMean}    = $refMean;
$ENV{refStdev}   = $refStdev;
$ENV{binWidth}   = $binWidth;
$ENV{maxBinSize} = $maxBinSize;
$ENV{minBinSize} = $minBinSize;
system("Rscript $ENV{LIB_DIR}/set_CN_emissions.R") == 0
    or die "error setting CN emission probabilities\n";

# create and solve HMM for each chromosome
print STDERR "processing chromosomes\n";
foreach my $CHROM(sort {$a cmp $b} @CHROMS){
    
    if ($ENV{IGNORE_CHROMS} =~ m/$CHROM,/) {
        my $msg = "ignoring $CHROM\n";
        print STDERR $msg;
        next;
    } 
    
    # load the bins
    print STDERR "  $CHROM\n";
    my $binH = getBinsH($CHROM, \my %stats);
    loadBins($binH);

    # solve the modal copy number HMM
    solveCopyNumberHMM($CHROM);
    
    # solve the per-sample CNV HMM 
    foreach my $sampleI(0..$#SAMPLES){
        solveCnvHMM($CHROM, $sampleI, $stats{weights});
    }    

    # write the updated bins file
    foreach my $bin(@bins){
        print  join("\t", @$bin)."\n";
    }
}

# HMM subs
sub solveCopyNumberHMM { # solve for the modal copy number of all samples
    my ($CHROM) = @_;
    my $ajdLenCol  = $binFileCols{adjLen};
    my $copyNumCol = $binFileCols{copyNum};

    # convert bin data to observation indices, correlated to bin indices
    my $tmpFile = "$ENV{TMP_DIR}/$GROUP_NAME.$CHROM.CN.segment.data";
    open(my $tmpH, ">", $tmpFile) or die "could not open $tmpFile for writing: $!\n";
    foreach my $binI(0..$#bins){ 
        print $tmpH join("\t", getBinI($bins[$binI][$ajdLenCol]), $binI), "\n"; 
    }
    close $tmpH;
    
    # convert resulting CN state indices to copy number
    open(my $hmmH, "-|", "cat $tmpFile | segment -e $eProbFile -z 0.1 -p 0.95")
        or die "could not open CN segment stream: $!\n";
    while (my $line = <$hmmH>) {
        chomp $line;
        my ($CN, $binI) = split("\t", $line);
        $CN > $MAX_CN and $CN = $MAX_CN;
        $bins[$binI][$copyNumCol] = $CN;
    }
    close $hmmH;
    unlink $tmpFile;    
}

sub solveCnvHMM { # solve for the copy number change of each sample
    my ($CHROM, $sampleI, $weights) = @_;
    my $SAMPLE     = $SAMPLES[$sampleI];
    $ENV{weight}   = $$weights[$sampleI]; # sample level data sent to R  
    my $medCovCol  = $binFileCols{medCov};    
    my $copyNumCol = $binFileCols{copyNum};    
    my $rawCovCol  = $nFileLeadCol + $sampleI;
    my $dCNCol     = $nFileLeadCol + @SAMPLES*2 + $sampleI;
    
    print STDERR "    $SAMPLE\n";
    
    # set the tmp files
    $ENV{dataFile}  = "$ENV{TMP_DIR}/$GROUP_NAME.$CHROM.$SAMPLE.CNV.segment.data";
    $ENV{eProbFile} = "$ENV{TMP_DIR}/$GROUP_NAME.$CHROM.$SAMPLE.CNV.emiss.prob.txt";
    
    # extract bin-level data for R, correlated to bin indices
    open(my $tmpH, ">", $ENV{dataFile}) or die "could not open $ENV{dataFile} for writing: $!\n";
    foreach my $binI(0..$#bins){     
        print $tmpH join("\t", (map {$bins[$binI][$_]}($medCovCol, $copyNumCol, $rawCovCol)), $binI), "\n"; 
    }
    close $tmpH;
    
    # calculate the bin-specific emission probabilities
    system("Rscript $ENV{LIB_DIR}/CNV_HMM.R") == 0
        or die "error setting CNV emission probabilities\n";
        
    # convert resulting cnDelta state indices to copy number changes
    open(my $hmmH, "-|", "cat $ENV{eProbFile} | segment -z 0.95 -p 0.95 -w")
        or die "could not open CNV segment stream: $!\n";
    while (my $line = <$hmmH>) {
        chomp $line;
        my ($dCN, $binI) = split("\t", $line);
        $bins[$binI][$dCNCol] = $dCN - 2; 
    } 
    close $hmmH;        
    unlink $ENV{dataFile};
    unlink $ENV{eProbFile};
}

# bin size subs
sub getBinSize {
    my ($varBinSize) = @_;
    return int($varBinSize / $binWidth + 0.5) * $binWidth; 
}
sub getBinI {
    my ($varBinSize) = @_;
    my $sizeBin = getBinSize($varBinSize);
    if ($sizeBin < $minBinSize ) {
        return 0;
    } elsif($sizeBin > $maxBinSize){
        return $maxBinI;
    } else {
        return ($sizeBin - $minBinSize) / $binWidth + 1;
    }
}
