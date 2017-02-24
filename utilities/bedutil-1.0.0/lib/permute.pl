#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'permute' takes feature BED lines on STDIN as input and prints:
#     1) the BED itself
#     2) $N_ITERATIONS random permutations of the BED
#     3) if desired, a single control iteration with the features offset
# in a single stream to STDOUT.  The iteration number is appended as the
# last field (==-1 for the input BED features, ==-2 for the offset control).  
# The output stream is suitable for input into any algorithm that can append
# a score field to the permuted BED features.  
#----------------------------------------------------------------------
# Options are:
#     CHROM_FILE         file with chromosome information
#                        format = "CHROM\tSIZE", one chromosome per line
#                        REQUIRED
#     EXCLUDE_FILES      comma-delimited list of BED files with excluded regions, e.g. gaps
#                        regions are always excluded on both strands and never padded
#                        OPTIONAL [default: no excluded regions]
#     MAX_EXCLUDED_FRAC  maximum fraction of a permuted feature allowed to overlap excluded regions
#                        excluded fraction is determined from the unpadded feature span
#                        OPTIONAL [default: 0, i.e. no part of feature may be excluded]
#     PROBES_FILE        BED file of probe positions used to restrict placement of permuted features 
#                        typical usage is for input features discovered by microarray
#                        OPTIONAL [default: probes are not considered]
#     MAX_PROBE_GAP      maximum allowed gap in bp between adjacent probes
#                        gaps larger than this size will be excluded from permuted feature placement
#                        OPTIONAL [default: 100000, i.e. 100 Kb]
#     MIN_PROBES         minimum number of probes required to be present in a permuted feature
#                        probe count is determined from the unpadded feature span
#                        OPTIONAL [default: 5]
#     MAX_OVERLAP_FRAC   maximum fraction of a permuted feature allowed to overlap other features
#                        overlap fraction is determined from the unpadded feature span
#                        OPTIONAL [default: 0, i.e. no feature overlap is allowed]
#     IS_FEMALE          boolean indicating whether the target genome is female
#                        if female, features are not allowed to be placed on chrY
#                        if male, chrX and chrY are weighted by half
#                        OPTIONAL [default: FALSE, male is assumed]
#     USE_MITO           boolean indicating whether to include the mitochondial chrM
#                        OPTIONAL [default: FALSE]  
#     STRANDED           boolean indicating whether (permuted) features are strand-specific
#                        if stranded, permuted features will be placed randomly onto a strand
#                        otherwise, any input strand field is retained in the permuted features
#                        OPTIONAL [default: FALSE]  
#     COLLAPSE           boolean indicating whether to collapse (permuted) features before printing
#                        when counting, the count will be in the penultimate column
#                        OPTIONAL [default: FALSE] 
#     N_ITERATIONS       the number of permutations to perform
#                        OPTIONAL [default: 100]
#     OFFSET_CONTROL     the number of bases by which to offset the control iteration
#                        all features are systematically offset by adding $OFFSET_CONTROL
#                        to their start and end coordinates
#                        OPTIONAL [default: 0, no offset control is created]
#----------------------------------------------------------------------
# The following options are inherited from collapse:
#     PADDING            bp added to each end of each feature prior to collapse.
#                        PADDING > 0 will thus merge adjacent but non-overlapping
#                        features into a single contiguous feature that includes
#                        the padded empty gap space between them.
#                        OPTIONAL [default: 0, i.e unpadded]  
#     UNPAD              boolean whether to remove padding from the flanks of collapsed regions
#                        OPTIONAL [default: 0, i.e collapsed features retain outside padding]  
#     COUNT              boolean indicating whether to count the number of collapsed features
#                        the count is appended as the last output column
#                        OPTIONAL [default: FALSE] 
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil permute CHROM_FILE=/path/to/chrom.bed IS_FEMALE
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 72, 66);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    CHROM_FILE
    EXCLUDE_FILES
    MAX_EXCLUDED_FRAC
    PROBES_FILE
    MAX_PROBE_GAP
    MIN_PROBES
    MAX_OVERLAP_FRAC
    IS_FEMALE
    USE_MITO 
    STRANDED 
    COLLAPSE
    N_ITERATIONS    
    OFFSET_CONTROL
    PADDING
    UNPAD
    COUNT
);
%booleanOptions = map { $_ => 1} qw (
    IS_FEMALE
    USE_MITO 
    STRANDED 
    COLLAPSE
    UNPAD
    COUNT
);
%defaultValues = (
    MAX_EXCLUDED_FRAC => 0,
    MAX_PROBE_GAP => 100000,
    MIN_PROBES => 5,
    MAX_OVERLAP_FRAC => 0,
    N_ITERATIONS => 100,
    OFFSET_CONTROL => 0
);
@requiredOptions = qw (
    CHROM_FILE
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
my $childHold = $ENV{IS_CHILD};
$ENV{IS_CHILD} = 1;
require "$scriptDir/collapse.pl";
$ENV{IS_CHILD} = $childHold;
#----------------------------------------------------------------------
my $reportStepSize = $ENV{N_ITERATIONS} / 10;
#######################################################################


#######################################################################
# collect the chromosome information
#----------------------------------------------------------------------
my ($n, $nChroms, $genomeSize, %chroms, %chromSizes) = (0, 0, 0);
open my $chromH, "<", $ENV{CHROM_FILE} or die "$error: could not open $ENV{CHROM_FILE}: $!\n";
while (my $line = <$chromH>){
    $line =~ m|^\s*#| and next;  # ignore comment lines
    chomp $line;
    $line =~ s/\r//g;
    my ($chrom, $chromSize) = split("\t", $line);
    $chrom or next;  # ignore blank lines
    my $ucChrom = uc($chrom);
    $ENV{IS_FEMALE} and $ucChrom eq "CHRY" and next;   
    $ENV{USE_MITO} or ($ucChrom eq "CHRM" and next);     
    my $error = "$error: invalid chromosome:\n$line\n";  
    parseInt(\$chromSize, $error);
    $nChroms++; 
    $chromSizes{$chrom} = $chromSize;
    my $nCopies;
    if($ENV{IS_FEMALE}){
        $nCopies = 1;
    } elsif ($ucChrom eq "CHRX" or $ucChrom eq "CHRY") {
        $nCopies = 1;
    } else {
        $nCopies = 2;
    }
    $genomeSize += ($nCopies * $chromSize);
    $chroms{$n++} = [$chrom, $chromSize];
    $nCopies == 2 and $chroms{$n++} = [$chrom, $chromSize];
}
close $chromH;
my $maxChromN = $n - 1;
$nChroms or die "$error: no chromosomes found in $ENV{CHROM_FILE}\n";
$genomeSize or die "$error: effective genome size was 0 for $ENV{CHROM_FILE}\n";
print STDERR "$feedback: $nChroms chromosomes extracted from $ENV{CHROM_FILE}\n";
print STDERR "$feedback: $genomeSize effective genome size\n";
#######################################################################


#######################################################################
# collect the excluded regions
#----------------------------------------------------------------------
my %tmp;
my @envHold = ($ENV{PADDING}, $ENV{STRANDED}, $ENV{COUNT});
($ENV{PADDING}, $ENV{STRANDED}, $ENV{COUNT}) = (0, 0, 0);  # excluded regions are never padded or strand specific
if($ENV{EXCLUDE_FILES}){
    my $nExcludedRegions = 0;
    foreach my $excludeFile(split(",", $ENV{EXCLUDE_FILES})){
        $excludeFile or next;
        loadBedFile($excludeFile, $error, \my$n, undef, \%tmp);
        $nExcludedRegions += $n;
    }
    print STDERR "$feedback: $nExcludedRegions excluded regions extracted from $ENV{EXCLUDE_FILES}\n"; 
}
#######################################################################   
    
    
#######################################################################
# collect the microarray probes
#----------------------------------------------------------------------
my (%probes, %probeCounts);
my $posIndexScalar = 100000;
if($ENV{PROBES_FILE}){
    loadBedFile($ENV{PROBES_FILE}, $error, \my$nProbes, \my@probes);
    $nProbes or die "$error: no probes found in $ENV{PROBES_FILE}\n";   
    print STDERR "$feedback: $nProbes probes extracted from $ENV{PROBES_FILE}\n";             
    foreach my $probe(@probes){ 
        my $position = $$probe[1];
        my $posIndex = int($position / $posIndexScalar);  # use a position index to improve speed with which probe counts are calculated
        push @{$probes{$$probe[0]}{$posIndex}}, $position;
    }
    my $nProbeGaps = 0;
    foreach my $chrom(keys %probes){  # exclude large gaps in probe coverage
        my $prevPosition = 0;    
        foreach my $posIndex(sort {$a <=> $b} keys %{$probes{$chrom}}){
            $probeCounts{$chrom}{$posIndex} = scalar(@{$probes{$chrom}{$posIndex}});
            foreach my $position(sort {$a <=> $b} @{$probes{$chrom}{$posIndex}}){ 
                $position - $prevPosition > $ENV{MAX_PROBE_GAP} and exludeProbeGap(\$nProbeGaps, $chrom, $prevPosition, $position + 1); 
                $prevPosition = $position;
            }
        }
        $chromSizes{$chrom} and $chromSizes{$chrom} - $prevPosition > $ENV{MAX_PROBE_GAP} and exludeProbeGap(\$nProbeGaps, $chrom, $prevPosition, $chromSizes{$chrom}); 
    }
    print STDERR "$feedback: $nProbeGaps gaps present in probe set\n"; 
}
sub exludeProbeGap {
    my ($nProbeGaps, $chrom, $start, $end) = @_;
    $$nProbeGaps++;
    push @{$tmp{$chrom}{0}{$start}{$end}}, [];
}
#######################################################################   
    
    
#######################################################################
# collapse the accumulated excluded regions
#----------------------------------------------------------------------
my ($excludedBases, %exclusions) = (0);
foreach my $chrom(keys %tmp){    
    collapseChromStrand($chrom, 0, \%tmp, \%exclusions);
    foreach my $exclusion(@{$exclusions{$chrom}{0}}){
        my ($start, $end) = @$exclusion;
        $excludedBases += ($end - $start);
    }
}     
print STDERR "$feedback: $excludedBases bases excluded from genome\n"; 
%tmp = ();
($ENV{PADDING}, $ENV{STRANDED}, $ENV{COUNT}) = @envHold;
#######################################################################


#######################################################################
# collect, validate and print the actual features
#----------------------------------------------------------------------
my $paddingHold = $ENV{PADDING};
$ENV{PADDING} = 0;  # load features without padding for assessing permutation placements (padding only applies to collapse)
loadBedFile(undef, $error, \our$nFeatures, \our@features, \our%features);
$ENV{PADDING} = $paddingHold;
$nFeatures or exit;
print STDERR "$feedback: $nFeatures features provided on STDIN\n";
$ENV{IS_CHILD} and return 1;
$ENV{OFFSET_CONTROL} and printPermutation(*STDOUT, offsetFeaturesHash(\%features), -2);
printPermutation(*STDOUT, \%features, -1);
runPermutation(*STDOUT);
#----------------------------------------------------------------------
sub printPermutation {
    my ($outH, $features, $j) = @_;
    if($ENV{COLLAPSE}){
        my $paddedFeatures = padFeaturesHash($features);
        collapseBED($paddedFeatures, \my%collapsed);
        printRegionsHash($outH, \%collapsed, undef, $j);  
    } else {
        printFeaturesHash($outH, $features, $j);   
    }
}
#######################################################################


#######################################################################
# randomly permute the input features
#----------------------------------------------------------------------
sub runPermutation {
    my ($outH, $suppressFeedback) = @_;     
    foreach my $j(0..($ENV{N_ITERATIONS} - 1)){
        my (%permFeatures, %permRegions);
        FEATURE: foreach my $feature(@features){
            my @feature = @$feature;
            my $featureSize = $feature[2] - $feature[1];
            REPEAT_FEATURE:
            my $baseIndex = int(rand $genomeSize);  # randomly place the feature
            my $cumChromSize = 0;
            foreach my $n(0..$maxChromN){
                my ($chrom, $chromSize) = @{$chroms{$n}};
                if($baseIndex < $cumChromSize + $chromSize){              
                    my $permChrom = $chrom;          
                    my $permStart = $baseIndex - $cumChromSize;        
                    my $permEnd = $permStart + $featureSize;  # end coordinate 1-indexed for BED format        
                    $permEnd > $chromSize and goto REPEAT_FEATURE;  # features cannot wrap chromosome ends 
                    isExcluded($permChrom, 0, $permStart, $permEnd, \%exclusions, $ENV{MAX_EXCLUDED_FRAC}, undef) 
                        and goto REPEAT_FEATURE;                    
                    my $permStrand = $ENV{STRANDED} ? int(rand 2) + 1 : 0;
                    my $permStrandIndex = $permStrand ?  ($permStrand == 1 ? '+' : '-') : $permStrand;                        
                    isExcluded($permChrom, $permStrandIndex, $permStart, $permEnd, \%permRegions, $ENV{MAX_OVERLAP_FRAC}, 1) 
                        and goto REPEAT_FEATURE;
                    hasSufficientProbes($permChrom, $permStart, $permEnd) or goto REPEAT_FEATURE;  
                    $ENV{STRANDED} and $feature[5] = $permStrandIndex; 
                    push @{$permFeatures{$permChrom}{$permStrandIndex}{$permStart}{$permEnd}}, [@feature[3..$#feature]];    
                    delete $permRegions{$permChrom}{$permStrandIndex};
                    collapseChromStrand($permChrom, $permStrandIndex, \%permFeatures, \%permRegions);
                    next FEATURE;
                }
                $cumChromSize += $chromSize; 
            }   
        }
        printPermutation($outH, \%permFeatures, $j);  
        $j and ($suppressFeedback or ($j % $reportStepSize or print STDERR "$j permutation iterations completed\n"));
    }
}
#----------------------------------------------------------------------
sub isExcluded {  # score permuted feature for excessive overlap with excluded regions
    my ($permChrom, $permStrandIndex, $permStart, $permEnd, $exclusions_, $maxExcludedFrac, $bidirectional) = @_;
    $$exclusions_{$permChrom} or return;
    my $exclusions = $$exclusions_{$permChrom}{$permStrandIndex};
    $exclusions or return;    
    $maxExcludedFrac >= 1 and return;
    my $permOverlap = 0;
    foreach my $exclusion(@$exclusions){    
        my ($exclusionStart, $exclusionEnd) = @$exclusion; 
        $exclusionStart < $permEnd and $exclusionEnd > $permStart or next;
        my $exclusionOverlap;        
        if($exclusionStart <= $permStart and $exclusionEnd >= $permEnd){
            $exclusionOverlap = $permEnd - $permStart;
        } elsif($exclusionStart >= $permStart and $exclusionEnd <= $permEnd) {
            $exclusionOverlap = $exclusionEnd - $exclusionStart;
        } elsif($exclusionStart < $permStart) {
            $exclusionOverlap = $exclusionEnd - $permStart;
        } elsif($exclusionEnd > $permEnd) {
            $exclusionOverlap = $permEnd - $exclusionStart; 
        } else {
            die "$error: isExcluded failed to match an overlap state\n";  # code error, should never occur
        }
        $bidirectional and $exclusionOverlap/($exclusionEnd - $exclusionStart) > $maxExcludedFrac and return 1;   
        $permOverlap += $exclusionOverlap;
    }
    return ($permOverlap/($permEnd - $permStart) > $maxExcludedFrac);
}
#----------------------------------------------------------------------
sub hasSufficientProbes {  # enforce minimum probe count
    my ($permChrom, $permStart, $permEnd) = @_;
    $ENV{PROBES_FILE} or return 1;
    my $probes = $probes{$permChrom};
    $probes or return undef; 
    $permEnd = $permEnd - 1;
    my $startIndex = int($permStart / $posIndexScalar); 
    my $endIndex = int($permEnd / $posIndexScalar); 
    my $N = 0;  
    if($endIndex - $startIndex > 1) { 
        foreach my $posIndex(($startIndex+1)..($endIndex-1)){
            my $n = $probeCounts{$permChrom}{$posIndex};
            $n and $N += $n;
            $N >= $ENV{MIN_PROBES} and return 1;
        }
        getJunctionCounts(\$N, $probes, $startIndex, $endIndex, $permStart, $permEnd) and return 1;
    } elsif($startIndex == $endIndex){
        getProbeCount(\$N, $probes, $startIndex, $permStart, $permEnd) and return 1;
    } else {
        getJunctionCounts(\$N, $probes, $startIndex, $endIndex, $permStart, $permEnd) and return 1;
    }
    return undef;
}
sub getJunctionCounts {
    my ($N, $probes, $startIndex, $endIndex, $permStart, $permEnd) = @_;
    getProbeCount($N, $probes, $startIndex, $permStart, ($startIndex + 1) * $posIndexScalar - 1) and return 1;  
    getProbeCount($N, $probes, $endIndex, $endIndex * $posIndexScalar, $permEnd) and return 1;  
}
sub getProbeCount {
    my ($N, $probes, $posIndex, $start, $end) = @_;
    foreach my $position(@{$$probes{$posIndex}}){
         $position >= $start or next;
         $position <= $end or next;
         $$N++;
         $$N >= $ENV{MIN_PROBES} and return 1;
    }
}
#######################################################################

1;

