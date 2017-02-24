use strict;
use warnings;

# find equally-weighted variable-width bins

# require common subs
require "$ENV{LIB_DIR}/common.pl";
use vars qw(@Ns %breaks @bins %binCols $nLeadCol);

# global variables
our $GROUP_NAME = $ENV{GROUP_NAME};
our @SAMPLES = split(" ", $ENV{SAMPLES});
our ($CHROM, $CHROM_SIZE, $STATS_FILE) = @ARGV;
if ($ENV{IGNORE_CHROMS} =~ m/$CHROM,/) {
    my $msg = "ignoring $CHROM\n";
    print $msg;
    print STDERR $msg;
    exit;
}

# recover the reference bin information
getStats($ENV{REF_STAT_FILE}, \our %refStats);
our $VAR_BIN_COUNT = $refStats{VAR_BIN_COUNT}[0];
our $REF_PLOIDY    = $refStats{REF_PLOIDY}[0];
print STDERR "using VAR_BIN_COUNT = $VAR_BIN_COUNT\n";

# parse to find bins
loadBreakMaps();
my $medianN = median(@Ns);
my @weights = map {$_ / $medianN} @Ns;
print STDERR "parsing bin ends\n";
setBinEnds();
%breaks = ();

# thread each sample to find its counts in each bin
our (@medians);
print STDERR "calculating bin counts for samples:\n";
foreach my $SAMPLE(@SAMPLES){
    print STDERR "    $SAMPLE\n";
    my $mapH = getMapH($SAMPLE, \my %stats);
    my $avgLen = $stats{avgLen}[0]; 
    our ($prevBreak, $cov, $binCount, $binI, @counts) = (0, 0, 0, 0);
    our $binStartCol = $binCols{start};    
    our $binEndCol   = $binCols{end};
    my $buffer; 
    while (read($mapH, $buffer, 4)) {
        my $break = unpack('L', $buffer);
        processBreak($break);
        $prevBreak = $break;
        read($mapH, $buffer, 2);        
        $cov += unpack('s', $buffer) / $avgLen;
        $binI > $#bins and last;
    }
    close $mapH;
    my $lastBreak = $CHROM_SIZE + 1;
    $lastBreak > $prevBreak and processBreak($lastBreak);
    sub processBreak {
        my ($break) = @_;    
        my $spanCount = ($break - $prevBreak) * $cov;
        $binCount += $spanCount; 
        while ($bins[$binI] and $break > $bins[$binI][$binEndCol] and $break > $prevBreak) {
            my $excess = $spanCount * ($break - $bins[$binI][$binEndCol] - 1) / ($break - $prevBreak);
            my $count = roundCount($binCount - $excess);
            push @{$bins[$binI]}, $count;
            push @counts, $count;
            $binCount = $excess;
            $binI++;
            $binI > $#bins and last;
            $prevBreak = $bins[$binI][$binStartCol];
            $spanCount = $binCount;   
        }   
    }
    push @medians, median(@counts);   
}

# calculated the adjusted bin sizes
print STDERR "calculating bin adjustments\n";
my $maxSampleI  = $#SAMPLES;
my $medCovCol   = $binCols{medCov};
my $adjLenCol   = $binCols{adjLen};
my $noGapLenCol = $binCols{noGapLen};
my $nBad = 0;
foreach my $bin(@bins){
    unless(defined $$bin[$nLeadCol]){ $nBad++; next }
    # normalize bin counts for sample depth
    my @normBinCounts =
        map {roundCount($$bin[$_ + $nLeadCol] / $weights[$_])} 0..$maxSampleI;
    push @$bin, @normBinCounts;         
    # the most typical sample count for this bin
    $$bin[$medCovCol] = roundCount(median(@normBinCounts)) || 1;
    # adjust bin size in case a CNV is in one sample
    $$bin[$adjLenCol] = int($$bin[$noGapLenCol] / ($$bin[$medCovCol] / $VAR_BIN_COUNT) + 0.5); 
}
print STDERR "    $nBad bad bins (with no data)\n";
print STDERR "calculating adjusted bin size distribution\n";
my ($binMedian, $binMean, $binStdev) = getBinStats();

# commit the stats
printStats($STATS_FILE,
    ['samples',      @SAMPLES],
    ['Ns',           @Ns], 
    ['weights',      @weights],
    ['medians',      @medians],
    ['adjBinMedian', $binMedian],
    ['adjBinMean',   $binMean], 
    ['adjBinStdev',  $binStdev] 
);

# commit the table of binned values for all samples
foreach my $bin(@bins){
    print join("\t", $CHROM, @$bin)."\n";
}
