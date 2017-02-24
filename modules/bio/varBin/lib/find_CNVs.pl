use strict;
use warnings;

# create a file optimized for plotting via MIBrowser

# require common subs
require "$ENV{LIB_DIR}/common.pl";
use vars qw(%plotCols $nPlotLeadCol);

# global variables
our $GROUP_NAME = $ENV{GROUP_NAME};
our @SAMPLES = split(" ", $ENV{SAMPLES});
my $MIN_CNV_BINS = $ENV{MIN_CNV_BINS};

# set input and output columns
my $chromCol = $plotCols{chrom};
my $startCol = $plotCols{start};
my $endCol   = $plotCols{end};
my $nSamples = @SAMPLES;
    
# thread data for runs of CNVs
print STDERR "finding runs of deviant bins by sample, i.e. CNVs\n";
my (%inRun, %bins);
while (my $line = <STDIN>) {
    chomp $line;
    my @bin = split("\t", $line);
    my $chrom = $bin[$chromCol];
    foreach my $i(0..$#SAMPLES){
        my $dcnHMMCol = $nPlotLeadCol + $nSamples*3 + $i;   
        my $cnv = $bin[$dcnHMMCol] < 0 ? 'loss' : ($bin[$dcnHMMCol] > 0 ? 'gain' : '');    
        if ($inRun{$i} and
            ($chrom ne $inRun{$i}[0] or
             $cnv   ne $inRun{$i}[1])){        
            commitCNV($i);
        }
        if ($cnv) {
            $inRun{$i} = [$chrom, $cnv];
            push @{$bins{$i}}, \@bin;
        }
    }
}
foreach my $i(0..$#SAMPLES){
    $inRun{$i} and commitCNV($i);
}

sub commitCNV {
    my ($i) = @_;
    my $cnvBins = $bins{$i};
    my $nBins = @$cnvBins; 
    if($nBins > $MIN_CNV_BINS) {
        my $chrom = $$cnvBins[0][$chromCol];
        my $start = $$cnvBins[0][$startCol];
        my $end   = $$cnvBins[$#{$cnvBins}][$endCol];
        my @smpDCN;
        foreach my $j(0..$#SAMPLES){
            my $dcnRawCol = $nPlotLeadCol + $nSamples*2 + $j;
            push @smpDCN, roundCount(mean(map {$$_[$dcnRawCol]} @$cnvBins), 100);
        }
        print join("\t", $chrom, $start, $end, "$chrom:$start-$end", $nBins, '.',
                         $inRun{$i}[1], $SAMPLES[$i], $i, @smpDCN), "\n";
    }
    delete $inRun{$i};
    delete $bins{$i};
}
