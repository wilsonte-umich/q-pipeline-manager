use strict;
use warnings;

# create a file optimized for plotting via MIBrowser

# require common subs
require "$ENV{LIB_DIR}/common.pl";
use vars qw(@bins %binFileCols $nFileLeadCol %plotCols $nPlotLeadCol);

# global variables
our $GROUP_NAME = $ENV{GROUP_NAME};
our @SAMPLES = split(" ", $ENV{SAMPLES});

# get the reference chromosome information
print STDERR "loading reference statistics\n";
getStats($ENV{REF_STAT_FILE}, \our %refStats);
my $REF_PLOIDY  = $refStats{REF_PLOIDY}[0];
my $refMean     = $refStats{refMean}[0];

# set input and output columns
my $medCovIn  = $binFileCols{medCov};    
my $adjLenIn  = $binFileCols{adjLen};
my $cnHMMIn   = $binFileCols{copyNum};
my $medCovOut = $plotCols{medCov};    
my $adjLenOut = $plotCols{adjLen};
my $cnHMMOut  = $plotCols{cnHMM};
my $cnRawOut  = $plotCols{cnRaw};
my $nSamples = @SAMPLES;
    
# fill the plot values
print STDERR "parsing plot values\n";
while (my $line = <STDIN>) {
    chomp $line;
    my @bin = split("\t", $line);
    my @plot = @bin[0..5];
    $plot[$medCovOut] = $bin[$medCovIn];
    $plot[$adjLenOut] = $bin[$adjLenIn];
    $plot[$cnHMMOut]  = $bin[$cnHMMIn];
    $plot[$cnRawOut] = roundCount($refMean / $bin[$adjLenIn] * $REF_PLOIDY, 100);    
    foreach my $sampleI(0..$#SAMPLES){
        my $rawCovIn  = $nFileLeadCol + $sampleI;
        my $nrmCovIn  = $nFileLeadCol + $nSamples*1 + $sampleI;
        my $dcnHMMIn  = $nFileLeadCol + $nSamples*2 + $sampleI;
        my $rawCovOut = $nPlotLeadCol + $sampleI;
        my $nrmCovOut = $nPlotLeadCol + $nSamples*1 + $sampleI;
        my $dcnRawOut = $nPlotLeadCol + $nSamples*2 + $sampleI;
        my $dcnHMMOut = $nPlotLeadCol + $nSamples*3 + $sampleI;
        $plot[$rawCovOut] = roundCount($bin[$rawCovIn], 10);
        $plot[$nrmCovOut] = roundCount($bin[$nrmCovIn], 10);
        $plot[$dcnHMMOut] = $bin[$dcnHMMIn];
        $plot[$dcnRawOut] = 
            roundCount(($bin[$nrmCovIn] - $bin[$medCovIn]) / $bin[$medCovIn] * $plot[$cnRawOut], 100);
    }
    print join("\t", @plot)."\n";
}
