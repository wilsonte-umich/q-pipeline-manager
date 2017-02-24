use strict;
use warnings;

# subs required by more than one varBin script

# global variables
use vars qw(@SAMPLES $CHROM $CHROM_SIZE $VAR_BIN_COUNT $GROUP_NAME $REF_PLOIDY %refStats);
our (@Ns, %breaks, @bins);

# define the bin file format, extended BED
our %binCols = ( # as held in memory, i.e. no chrom
    chrom       => -1,      
    start       =>  0,
    end         =>  1, 
    name        =>  2, # unused
    length      =>  3,
    strand      =>  4, # unused
    noGapLen    =>  5, 
    medCov      =>  6,
    adjLen      =>  7,
    copyNum     =>  8,
    # samples raw counts
    # samples normalized counts
    # samples delta copy number HMM 
);
our %binFileCols = map { $_ => $binCols{$_} + 1 } keys %binCols;
our $nLeadCol = 9;
our $nFileLeadCol = $nLeadCol + 1;

# define the plot file format
our %plotCols = ( # as held in memory, i.e. no chrom
    chrom       =>  0,      
    start       =>  1,
    end         =>  2, 
    name        =>  3, # unused
    medCov      =>  4, 
    strand      =>  5, # unused
    adjLen      =>  6, 
    cnRaw       =>  7,
    cnHMM       =>  8,
    # samples raw counts
    # samples normalized counts
    # samples delta copy number
    # samples delta copy number HMM 
);
our $nPlotLeadCol = 9;

# recover all stats from a stats file
sub getStats {
    my ($statsFile, $stats) = @_;
    open(my $statH, "<", $statsFile) or die "could not open $statsFile for reading: $!\n";
    while (my $line = <$statH>) {
        chomp $line;
        my ($name, @vals) = split("\t", $line);
        $$stats{$name} = \@vals;
    }
    close $statH;
}
sub printStats {
    my ($statsFile, @stats) = @_;
    open(my $statH, ">", $statsFile) or die " could not open $statsFile for writing: $!\n";
    foreach my $stat(@stats){
        my ($name, @vals) = @$stat;
        my $line = join("\t", $name, @vals)."\n";
        print $statH $line;
        print STDERR $line;        
    }
    close $statH;
}

# get file handle into break coverage map binary file
sub getMapH {
    my ($SAMPLE, $stats) = @_;
    my $mapFile = "$ENV{MAPS_DIR}/$SAMPLE/$SAMPLE.$CHROM.map";
    open(my $mapH, "-|", "slurp $mapFile") or die "could not open $mapFile for reading: $!\n";
    $stats and getStats("$mapFile.stats", $stats);
    return $mapH;
}

# merge the coverage maps of all input samples
sub loadBreakMaps {
    print STDERR "loading break maps for all samples\n";
    %breaks = @Ns = ();
    foreach my $SAMPLE(@SAMPLES){
        my $mapH = getMapH($SAMPLE, \my %stats);
        push @Ns, $stats{nFrags}[0];
        my $avgLen = $stats{avgLen}[0];
        my $buffer;
        while (read($mapH, $buffer, 4)) {        
            my $break = unpack('L', $buffer);
            read($mapH, $buffer, 2);        
            $breaks{$break} # sum the coverage changes at all break points across all samples
                = ($breaks{$break} || 0)
                +  unpack('s', $buffer) / $avgLen;
        }
        close $mapH;
    }    
}

# thread the merged breaks to find the optimal bin endpoints
sub setBinEnds {
    @bins = ();
    my $SUM_BIN_COUNT = $VAR_BIN_COUNT * @SAMPLES;
    my ($prevBreak, $cov, $binCount, $binStart) = (0) x 10;
    BREAK: foreach my $break(sort {$a <=> $b} keys %breaks){
        my $spanCount = ($break - $prevBreak) * $cov;
        $binCount += $spanCount;
        while ($binCount >= $SUM_BIN_COUNT and $spanCount) {
            my $binEnd = int(0.5 + $prevBreak - 1 +
                             ($break - $prevBreak) * (1 - ($binCount - $SUM_BIN_COUNT) / $spanCount));
            push @bins, [$binStart, $binEnd, '.', $binEnd - $binStart + 1, '.',
                         'na', 'na', 'na', 'na'];
            $binStart = $binEnd + 1;        
            $binCount = ($break - $binStart) * $cov;
            $prevBreak = $binStart;
            $spanCount = $binCount;
        }
        $prevBreak = $break;
        $cov += $breaks{$break};   
    }
    if ($binStart <= $CHROM_SIZE) {
        push @bins, [$binStart, $CHROM_SIZE, '.', $CHROM_SIZE - $binStart + 1, '.',
                     'na', 'na', 'na', 'na'];
    }
    correctForGaps();
}

# correct bin lengths by substracting out contained gap lengths
sub correctForGaps {
    my $tmpFile = "$ENV{TMP_DIR}/$GROUP_NAME.$CHROM.gap.tmp.bed";
    open(my $outH, ">", $tmpFile) or die "could not open $tmpFile for writing: $!\n";
    my @cols = ($binCols{start}, $binCols{end});
    foreach my $bin(@bins){
        print $outH join("\t", $CHROM, @$bin[@cols]), "\n";
    }
    close $outH;
    open(my $inH, "-|",
        "bedtools intersect -wao -a $tmpFile -b $ENV{GAP_FILE} | ".
        "awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$NF}' | ".
        "groupBy -g 1,2,3 -c 4 -o sum | cut -f 4"
    ) or die "could not open correctForGaps stream: $!\n";
    my $i = 0;
    my $lenCol      = $binCols{length};
    my $noGapLenCol = $binCols{noGapLen};
    while (my $gapLen = <$inH>) {
        chomp $gapLen;
        $gapLen or $gapLen = 0;
        $bins[$i][$noGapLenCol] = $bins[$i][$lenCol] ? $bins[$i][$lenCol] - $gapLen : 0;
        $i++;
    }
    close $inH;
    unlink $tmpFile;
}

# return information on the main peak of bins
sub getRefStats {
    my ($start, $end) = @_;
    my @sizes;
    my $startCol = $binCols{start};
    my $endCol   = $binCols{end};
    my $lenCol   = $binCols{noGapLen};
    foreach my $bin(@bins){
        $$bin[$startCol] >= $start or next;
        $$bin[$endCol]   <= $end   or next;
        push @sizes, $$bin[$lenCol];
    }
    my $median = median(@sizes);
    my $min = $median - $median / 2;
    my $max = $median + $median / 2;
    my @peak;    
    foreach my $size(@sizes){
            $size >= $min
        and $size <= $max
        and push @peak, $size;
    }
    return (median(@peak), stdev(@peak)); # median, mean, stdev
}
sub getBinStats {
    my $refMean  = $refStats{refMean}[0];
    my $refStdev = $refStats{refStdev}[0];
    my $nSD = 2.5;
    my %sizes;
    my $lenCol   = $binCols{adjLen};
    foreach my $CN(1..$REF_PLOIDY + 1){
        my ($cnMean, $cnStdev) = adjustSizeStats($refMean, $refStdev, $CN);
        my $min = $refMean - $refStdev * $nSD;
        my $max = $refMean + $refStdev * $nSD;        
        foreach my $bin(@bins){
                $$bin[$lenCol] >= $min
            and $$bin[$lenCol] <= $max
            and push @{$sizes{$CN}}, $$bin[$lenCol];
        }
    }
    my @nullChrom = (0, 0, 0);
    my @availCNs = keys %sizes;
    @availCNs or return @nullChrom; 
    my @CNs = sort { @{$sizes{$b}} <=> @{$sizes{$a}} } @availCNs;    
    my $CN = $CNs[0];
    @{$sizes{$CN}} > 100 or return @nullChrom; 
    my @peak = @{$sizes{$CN}};
    return (median(@peak), stdev(@peak)); # median, mean, stdev
}

# get file handle into a chromosome's bin file
sub getBinsH {
    my ($CHROM, $stats) = @_;
    my $binsFile = "$ENV{BINS_DIR}/$GROUP_NAME.$CHROM.bins.bgz";
    open(my $binH, "-|", "slurp $binsFile | bgzip -cd")
        or die "could not open $binsFile for reading: $!\n";
    $stats and getStats("$binsFile.stats", $stats);
    return $binH;
}

# load all of the bins from a file (includes chrom column)
sub loadBins {
    my ($inH) = @_;
    @bins = ();
    while (my $line = <$inH>) {
        chomp $line;
        my @bin = split("\t", $line);
        defined $bin[$nFileLeadCol] or next;
        push @bins, \@bin;
    }
    close $inH;
}

# project the bin width distribution to a different copy number state
sub adjustSizeStats {
    my ($refMean, $refStdev, $CN) = @_;
    my $scalar = $REF_PLOIDY / $CN;
    return ($refMean * $scalar, $refStdev * $scalar); 
}

# math functions
sub median {
    my (@data) = sort {$a <=> $b} @_;
    my $i = @data / 2;
    my $upper = $data[$i];
    @data % 2 and return $upper;
    my $lower = $data[$i - 1];    
    return($lower + ($upper - $lower) / 2);
}
sub mean{
    my (@data) = @_;
    @data == 0 and die "no values sent to mean\n";
    my $sum = 0;
    foreach (@data) { $sum += $_ }
    return $sum / @data;
}
sub stdev{
    my (@data) = @_;
    @data == 1 and return 0;
    my $mean = mean(@data);
    my $sqSum = 0;
    foreach(@data) { $sqSum += ($mean-$_) ** 2 }
    return ($mean, ($sqSum / (@data-1)) ** 0.5);
}
sub roundCount {
    my ($count, $scalar) = @_;
    $scalar or $scalar = 1000;
    int($count * $scalar + 0.5) / $scalar;
}

1;
