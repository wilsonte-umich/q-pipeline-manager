use strict;
use warnings;

# optimize bin width for comparing copy number $REF_PLOIDY vs. $REF_PLOIDY + 1

# validated observations:
#   median, mean and stdev bin size are ~linear functions of $VAR_BIN_COUNT
#   but scale slightly differently so that CV is non-linear inverse to $VAR_BIN_COUNT
#   distribution of bin sizes in main peak is Normal and can be modeled as such

# require common subs
require "$ENV{LIB_DIR}/common.pl";
use vars qw(%binCols);

# global variables
our $GROUP_NAME = $ENV{GROUP_NAME};
our @SAMPLES = split(" ", $ENV{SAMPLES});
my $maxSamples = 10; # limit number of samples for time efficiency, assume all are similar
@SAMPLES > $maxSamples and @SAMPLES = @SAMPLES[0..$maxSamples-1];
our ($REF_REGION, $REF_STAT_FILE) = @ARGV;
our ($CHROM, $span) = split(':', $REF_REGION);
our ($START, $END)  = $span ? split('-', $span) : (1, 1e9);
our $REF_PLOIDY = $ENV{REF_PLOIDY};
our $CHROM_SIZE = 0;

# load the coverage map for all samples
loadBreakMaps();

# determine the optimum nominal bin count
our $VAR_BIN_COUNT;
my $gainScalar = $REF_PLOIDY / ($REF_PLOIDY + 1); # test counts for 1N gain from REF_PLOIDY
my $minNSD = 1.25; # minimum number of standard deviations to ~crossing point from each mean
my ($vbcInc, $vbc) = (25, 0);
do { $vbc += $vbcInc } until checkVBC();
sub checkVBC { # mimics what an average SINGLE sample analysis would look like
    $VAR_BIN_COUNT = int($vbc / @SAMPLES + 0.5);
    setBinEnds();
    my ($median, $mean, $stdev) = getRefStats($START, $END);
    my $meanGain  = $mean  * $gainScalar;
    my $stdevGain = $stdev * $gainScalar;
    print STDERR join("\t", $vbc, int($mean), int($stdev),
                                  int($meanGain), int($stdevGain)), "\n";
    return ($mean - $stdev * $minNSD) > ($meanGain + $stdevGain * $minNSD);
}

# commit the selected value
$VAR_BIN_COUNT = $vbc;
setBinEnds();
my ($refMedian, $refMean, $refStdev) = getRefStats($START, $END);
printStats($REF_STAT_FILE,
    ['REF_REGION',    $REF_REGION],
    ['REF_PLOIDY',    $REF_PLOIDY], 
    ['VAR_BIN_COUNT', $VAR_BIN_COUNT],
    ['refMedian',     $refMedian],
    ['refMean',       $refMean], 
    ['refStdev',      $refStdev] 
);
