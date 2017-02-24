#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'crossing' takes feature BED lines on STDIN and calculates 
#  a score against a provided reference BED file based on whether the 
# input feature crosses one or more reference features.  Different score
# type options are summarized below. The BED lines are printed to STDOUT
# with the score appended as the last column.
#----------------------------------------------------------------------
# Options are:
#     REF_FILE           BED file of reference features to which input features are compared
#                        REQUIRED
#     SCORE_TYPE         the type of crossing score reported for each feature
#                        recognized score types are summarized below
#                        OPTIONAL [default: COUNT] 
#     MIN_CROSSING_FRAC  minimum crossing fraction required to be scored as boolean positive
#                        OPTIONAL [default: 1/input length, i.e. 1 bp overlap] 
#     PADDING            bp added to each end of each input and reference feature prior to comparison
#                        PADDING > 0 will thus find crossings between nearby adjacent features
#                        fractional scores are reported as a fraction of the padded input length
#                        but the features themselves are always unpadded in the output
#                        OPTIONAL [default: 0, i.e unpadded] 
#     STRANDED           boolean indicating whether features are strand-specific
#                        if stranded, only features on the same strand are considered matching
#                        OPTIONAL [default: FALSE] 
#     MAX_NULL           numeric score to report when the input feature matched no reference features
#                        applies to MAX_REF_LENGTH, MAX_REF_SCORE, WEIGHTED2, and DISTANCE_<TYPE>
#                        OPTIONAL [default: 'NULL', i.e. no value to report]
#----------------------------------------------------------------------
# Score types are:
#     COUNT              the number of reference features crossed by the input feature
#     COUNT_ENDS         the number of reference features into which one or both input features ends fell
#     LENGTH             the number of input bases that crossed reference feature(s)
#     FRACTION           the fraction of the input feature that crossed reference feature(s)
#     BOOLEAN            whether or not FRACTION exceeded MIN_CROSSING_FRAC (0,1)
#     BOOLEAN_ENDS       whether or not one or both input features ends fell inside a reference feature
#     WEIGHTED           FRACTION multiplied by the reference feature score field (column 5)
#                        summed over all reference features that were crossed
#                        implies score=zero for input portions not crossing a reference feature
#     WEIGHTED2          like WEIGHTED, except that FRACTION is calculated using only those
#                        portions of the input feature that actually crossed reference features
#                        if no reference features were crossed, the result is MAX_NULL                      
#     MAX_REF_LENGTH     the maximum length of all reference features crossed
#     MAX_REF_SCORE      the maximum score (BED column 5) of all reference features crossed
#     DISTANCE_CENTER    the distance in bp from the input feature to the closest reference feature
#     DISTANCE_EDGE      DISTANCE_CENTER is measured between centers, DISTANCE_EDGE between closest edges
#                        DISTANCE_EDGE = 0 for overlapping input and reference features
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil crossing REF_FILE=/path/to/ref.bed STRANDED
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 55, 49);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    REF_FILE
    SCORE_TYPE
    MIN_CROSSING_FRAC
    PADDING    
    STRANDED 
    MAX_NULL
);
%booleanOptions = map { $_ => 1} qw (
    STRANDED
);
%defaultValues = (
    SCORE_TYPE => 'COUNT',
    PADDING => 0,
    MAX_NULL => 'NULL'
);
@requiredOptions = qw (
    REF_FILE
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
my $childHold = $ENV{IS_CHILD};
$ENV{IS_CHILD} = 1;
require "$scriptDir/collapse.pl";
$ENV{IS_CHILD} = $childHold;
#----------------------------------------------------------------------
my @nonNullableScores = qw(COUNT COUNT_ENDS LENGTH FRACTION BOOLEAN BOOLEAN_ENDS WEIGHTED);
my @nullableScores =    qw(MAX_REF_LENGTH MAX_REF_SCORE WEIGHTED2 DISTANCE_CENTER DISTANCE_EDGE);
my %defaultScores;
foreach my $scoreType(@nonNullableScores){ $defaultScores{$scoreType} = 0 }
foreach my $scoreType(@nullableScores){    $defaultScores{$scoreType} = $ENV{MAX_NULL} }
my %initialScores = %defaultScores;
$initialScores{WEIGHTED2} = 0;
defined $defaultScores{$ENV{SCORE_TYPE}} or die "$error: unrecognized score type: $ENV{SCORE_TYPE}\n";
#######################################################################


#######################################################################
# collect and process the reference features
#----------------------------------------------------------------------
my (%refRegions);
loadBedFile($ENV{REF_FILE}, $error, \my$nRefFeatures, undef, \my%tmp);  # will pad the reference features
my @envHold = ($ENV{PADDING}, $ENV{COUNT});
($ENV{PADDING}, $ENV{COUNT}) = (0, 0);  # no need to pad when collapsing the reference features   
my ($collRefBases, $maxCollRefLength, $maxCollRefEnd, $nCollRefFeatures) = (0, 0, 0, 0);
collapseRefFeatures();
my $aveCollRefLength = int($collRefBases / $nCollRefFeatures);
my $refBinSize = int($maxCollRefEnd / 100);  # nominally split largest chromosome into 100 bins for faster searching (but higher memory use)
my $minRefBinSize = $aveCollRefLength * 2;
$refBinSize = $refBinSize > $minRefBinSize ?  $refBinSize : $minRefBinSize;   # but use larger bins if reference features are especially long
binRefFeatures();
print STDERR "$feedback: $nRefFeatures reference features loaded from $ENV{REF_FILE}\n"; 
print STDERR "$feedback: $nCollRefFeatures reference regions after collapsing\n"; 
print STDERR "$feedback: $collRefBases bases present in reference regions, including padding\n"; 
print STDERR "$feedback: longest reference region: $maxCollRefLength bases\n";
print STDERR "$feedback: average reference region: $aveCollRefLength bases\n";
($ENV{PADDING}, $ENV{COUNT}) = @envHold;
#######################################################################


#######################################################################
# collect BED features from STDIN and score
#----------------------------------------------------------------------
$ENV{IS_CHILD} and return 1;
loadBedFile(undef, $error, \my$nFeatures, undef, \my%features);
$nFeatures or exit;
print STDERR "$feedback: $nFeatures query features provided on STDIN\n";
scoreCrossing(*STDOUT, \%features, \%refRegions);
#######################################################################


#######################################################################
# worker subs
#----------------------------------------------------------------------
sub collapseRefFeatures {  # convert reference features to region hash, collapsing for certain score types
    my %collapse = map { $_ => 1 } qw (LENGTH FRACTION BOOLEAN BOOLEAN_ENDS);
    foreach my $chrom(keys %tmp){    
        foreach my $strandIndex(keys %{$tmp{$chrom}}){ 
            if($collapse{$ENV{SCORE_TYPE}}){ # collapse reference features to facilitate length, fraction and derived boolean scoring                
                collapseChromStrand($chrom, $strandIndex, \%tmp, \%refRegions);  
            } else {
                featureHashToRegionHash($chrom, $strandIndex, \%tmp, \%refRegions);  # do not collapse when counting, weighting, or taking max value       
            }
            $refRegions{$chrom}{$strandIndex} or next;   
            foreach my $refRegion(@{$refRegions{$chrom}{$strandIndex}}){  # collect reference feature stats
                my ($start, $end) = @$refRegion;
                my $featureLength = $end - $start;
                $collRefBases += $featureLength;
                $maxCollRefLength >= $featureLength or $maxCollRefLength = $featureLength;
                $maxCollRefEnd >= $end or $maxCollRefEnd = $end;
                $nCollRefFeatures++;
            }
        }
    }     
}
sub binRefFeatures {  # split reference features into bins on each chromosome strand to speed up crossing search
    $ENV{SCORE_TYPE} =~ m/^DISTANCE_/ and return;  # cannot bin distance assessments (at least not DISTANCE_CENTER), must check all on chrom/strand
    my %binned;
    foreach my $chrom(keys %refRegions){   
        foreach my $strandIndex(keys %{$refRegions{$chrom}}){
            foreach my $refRegion(@{$refRegions{$chrom}{$strandIndex}}){
                my ($start, $end) = @$refRegion;
                my ($startBin, $endBin) = getCrossedBins($start, $end);
                for (my $bin = $startBin; $bin <= $endBin; $bin += $refBinSize){
                    push @{$binned{$chrom}{$strandIndex}{$bin}}, $refRegion;  # reference feature is stored in every bin it crosses
                }
            }
        }
    } 
    %refRegions = %binned;
}
sub getCrossedBins {  # the lowest and highest bins crossed by a region
    my ($start, $end) = @_;
    return (getBinIndex($start), getBinIndex($end - 1));  # end bins are converted to 0-referenced start coordinates
}
sub getBinIndex {  # bin index is the lowest 0-referenced base number in the bin
    my ($position) = @_;
    return int($position / $refBinSize) * $refBinSize;
}
#----------------------------------------------------------------------
sub scoreCrossing {  # iterate over all input features and print with appended score
    my ($outH, $test, $reference) = @_;
    foreach my $chrom(keys %$test){   
        foreach my $strandIndex(keys %{$$test{$chrom}}){
            foreach my $start(keys %{$$test{$chrom}{$strandIndex}}){ 
                foreach my $end(keys %{$$test{$chrom}{$strandIndex}{$start}}){              
                    my $score = getCrossingScore($chrom, $strandIndex, $start, $end, $reference);
                    defined $score or $score = $defaultScores{$ENV{SCORE_TYPE}};  
                    foreach my $trailing(@{$$test{$chrom}{$strandIndex}{$start}{$end}}){ 
                        my @feature = ($chrom, $start, $end, @$trailing, $score);
                        unpadFeature(\@feature, 1);
                        print $outH join("\t", @feature), "\n"; 
                    } 
                }  
            }
        }
    }  
}
sub getCrossingScore {  # tabulate features overlap(s) with reference region(s); this sub is where scores get defined
    my ($chrom, $strandIndex, $testStart, $testEnd, $reference) = @_;
    $$reference{$chrom} or return;
    my $refRegions = $$reference{$chrom}{$strandIndex};
    $refRegions or return;
    my $testLength = $testEnd - $testStart;
    $testLength or return;    
    $ENV{SCORE_TYPE} =~ m/^DISTANCE_/ and return getMinDistance($testStart, $testEnd, $testLength, $refRegions);
    my ($startBin, $endBin) = getCrossedBins($testStart, $testEnd); 
    my %score = %initialScores;       
    my %encountered;
    for (my $bin = $startBin; $bin <= $endBin; $bin += $refBinSize){  # check all bins crossed by query feature
        my $binRegions = $$refRegions{$bin};
        $binRegions or next;
        foreach my $refRegion(@$binRegions){  
            $encountered{$refRegion} and next;  # only score a reference region once, even if it crossed multiple query bins
            $encountered{$refRegion}++;
            my ($refStart, $refEnd, $refName, $refScore) = @$refRegion; 
            $refStart <= $testEnd and $refEnd >= $testStart or next;
            $score{COUNT}++;
            (($testStart >= $refStart and $testStart <= $refEnd) or
             ($testEnd   >= $refStart and $testEnd   <= $refEnd)) and $score{COUNT_ENDS}++;
            my $refLength = $refEnd - $refStart; 
            ($score{MAX_REF_LENGTH} and $score{MAX_REF_LENGTH} ne $ENV{MAX_NULL} and $score{MAX_REF_LENGTH} >= $refLength)
                or $score{MAX_REF_LENGTH} = $refLength;
            $refScore or $refScore = 0;               
            if($ENV{SCORE_TYPE} eq 'MAX_REF_SCORE'){
                ($score{MAX_REF_SCORE} and $score{MAX_REF_SCORE} ne $ENV{MAX_NULL} and $score{MAX_REF_SCORE} >= $refScore)
                    or $score{MAX_REF_SCORE} = $refScore;
            }
            if($refStart <= $testStart and $refEnd >= $testEnd){
                scoreCrossingOverlap(\%score, $testLength, $testLength, $refScore);
            } elsif($refStart >= $testStart and $refEnd <= $testEnd) {
                scoreCrossingOverlap(\%score, $refLength, $testLength, $refScore);
            } elsif($refStart < $testStart) {
                scoreCrossingOverlap(\%score, $refEnd - $testStart, $testLength, $refScore);
            } elsif($refEnd > $testEnd) {
                scoreCrossingOverlap(\%score, $testEnd - $refStart, $testLength, $refScore);
            } else {
                die "$error: getCrossingScore failed to match an overlap state\n";  # code error, should never occur
            }
        }
    }
    scalar(keys %encountered) or return;
    $score{FRACTION} and $score{BOOLEAN} = $score{FRACTION} >= ($ENV{MIN_CROSSING_FRAC} ? $ENV{MIN_CROSSING_FRAC} : 1 / $testLength);
    $score{BOOLEAN} or $score{BOOLEAN} = 0;
    $score{BOOLEAN_ENDS} = $score{COUNT_ENDS} ? 1 : 0;
    $ENV{SCORE_TYPE} eq 'WEIGHTED2' and $score{WEIGHTED2} = ($score{FRACTION} ? $score{WEIGHTED2} / $score{FRACTION} : undef);
    return $score{$ENV{SCORE_TYPE}};
}
sub scoreCrossingOverlap {
    my ($score, $overlap, $testLength, $refScore) = @_;
    $$score{LENGTH} += $overlap;
    my $fraction = $overlap / $testLength;    
    $$score{FRACTION} += $fraction;       
    $ENV{SCORE_TYPE} eq 'WEIGHTED'  and $$score{WEIGHTED}  += $fraction * $refScore;  
    $ENV{SCORE_TYPE} eq 'WEIGHTED2' and $$score{WEIGHTED2} += $fraction * $refScore;   
}
sub getMinDistance {
    my ($testStart, $testEnd, $testLength, $refRegions) = @_;
    my $minDistance = 1e9;
    if($ENV{SCORE_TYPE} eq 'DISTANCE_CENTER'){
        my $testCenter = int(($testStart + $testLength/2) + 0.5);
        foreach my $refRegion(@$refRegions){
            my ($refStart, $refEnd) = @$refRegion;
            my $refCenter = int(($refStart + ($refEnd - $refStart)/2) + 0.5);
            my $distance = abs($testCenter - $refCenter);
            $minDistance <= $distance or $minDistance = $distance;
            $minDistance == 0 and return 0;
        }
    }else{
        foreach my $refRegion(@$refRegions){
            my ($refStart, $refEnd) = @$refRegion;
            $refStart <= $testEnd and $refEnd >= $testStart and return 0;  # crossing features return DISTANCE_EDGE = 0
            my $distance = $refEnd > $testStart ? $refStart - $testEnd : $testStart - $refEnd;
            $minDistance <= $distance or $minDistance = $distance;
            $minDistance == 0 and return 0;
        }    
    }
    return $minDistance;
}
#######################################################################

1;

