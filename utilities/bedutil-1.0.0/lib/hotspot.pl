#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'hotspot' takes feature BED lines on STDIN and collapses them to 
# create cluster regions of overlapping features.  The input features are 
# further subjected to permutation, with the permuted features subjected 
# to collapse within each iteration.  Next, the input features that 
# correspond to the cluster(s) with the greatest number of features are
# removed from the data, the correponding clusters are exluded from
# placement of permuted features, and the permutation is run again. This
# process of removing the highest ranking clusters and repeating the 
# permutation is repeated until no features remain.  p-value estimates are 
# made for each permutation as the fraction of iterations that had at least 
# as many features in clusters with at least as many features as the highest 
# ranking cluster(s) as were contained in 1,2,3... of the highest ranking 
# cluster(s). The logic of this approach is that the highest ranking clusters 
# removed at each round are assumed to be non-random and therefore would  
# confound the assessment of clusters with smaller numbers of features since  
# the problem is hypergeometric, that is, population sampling from a finite 
# population without replacement.
#----------------------------------------------------------------------
# The following options are inherited from permute:
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
#     IS_FEMALE          boolean indicating whether the target genome is female
#                        if female, features are not allowed to be placed on chrY
#                        if male, chrX and chrY are weighted by half
#                        OPTIONAL [default: FALSE, male is assumed]
#     USE_MITO           boolean indicating whether to include the mitochondial chrM
#                        OPTIONAL [default: FALSE]  
#     STRANDED           boolean indicating whether (permuted) features are strand-specific
#                        if stranded, permuted features will be placed randomly onto a strand
#                        and only features on the same strand are collapsed together
#                        OPTIONAL [default: FALSE]  
#     N_ITERATIONS       the number of permutations to perform in each permutation cycle
#                        OPTIONAL [default: 100]
#----------------------------------------------------------------------
# The following options are inherited from collapse:
#     PADDING            bp added to each end of each feature prior to collapse.
#                        PADDING > 0 will thus merge adjacent but non-overlapping
#                        features into a single contiguous feature that includes
#                        the padded empty gap space between them.
#                        OPTIONAL [default: 0, i.e unpadded]  
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil hotspot CHROM_FILE=/path/to/chrom.bed IS_FEMALE
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 67, 61);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    CHROM_FILE
    EXCLUDE_FILES
    MAX_EXCLUDED_FRAC
    PROBES_FILE
    MAX_PROBE_GAP
    MIN_PROBES 
    IS_FEMALE
    USE_MITO
    STRANDED
    N_ITERATIONS
    PADDING    
);
%booleanOptions = map { $_ => 1} qw (
    IS_FEMALE
    USE_MITO
    STRANDED
);
%defaultValues = (
);
@requiredOptions = qw (
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
$ENV{MAX_OVERLAP_FRAC} = 1;  # the nature of hotspot finding allows features to overlap in clusters
$ENV{COLLAPSE} = 1;
$ENV{COUNT} = 1;
$ENV{UNPAD} = 1;
$ENV{IS_CHILD} and return 1;
$ENV{IS_CHILD} = 1;
require "$scriptDir/permute.pl";  # permute.pl will load the unpadded features
#######################################################################


#######################################################################
# collect and collapse BED features loaded by permute.pl
#----------------------------------------------------------------------
use vars qw(%exclusions $nFeatures @features %features);
open my $actualCollapsedH, ">", \my$actualCollapsed or die $!;
printPermutation($actualCollapsedH, \%features, 0);  # collapse.pl will apply paddinng only during collapse
close $actualCollapsedH;
loadClusters(\$actualCollapsed, \my$nActualClusters, \my@actualClusters);
print STDERR "$feedback: $nActualClusters clusters resulted from collapsing the input features\n";
#######################################################################


#######################################################################
# run the permutation cycles
#----------------------------------------------------------------------
print "\t\t\t\t".
      "frequency >=N matching clusters\n";
print "remaining features\tfeatures in highest ranking cluster(s)\t# clusters\t# features\t".
      "1\t2\t3\t...\n";   
my $nRemainingFeatures = $nFeatures;      
my @workingClusters = @actualClusters;
my ($maxNFeatures, $nMaxClusters, $nMaxFeatures) = getMaxNFeatures(\my@maxNClusters, \my@remainingClusters);
while($maxNFeatures){  # step down from the most significant clusters, exluding them as the process continues
    runHotspotPermutation();
    $nRemainingFeatures = excludeMaxNClusters();
    @workingClusters = @remainingClusters;
    ($maxNFeatures, $nMaxClusters, $nMaxFeatures) = getMaxNFeatures(\@maxNClusters, \@remainingClusters);
}
print "\n";
#######################################################################


#######################################################################
# worker subs
#----------------------------------------------------------------------
sub loadClusters {
    my ($collapsed, $nClusters, $clusters) = @_;
    my $paddingHold = $ENV{PADDING};
    $ENV{PADDING} = 0;
    loadBedFile($collapsed, $error, $nClusters, $clusters);
    $ENV{PADDING} = $paddingHold;
}
#----------------------------------------------------------------------
sub getMaxNFeatures {  # find the highest remaining ranking clusters, i.e. with the greatest number of contained features
    my ($maxNClusters, $remainingClusters) = @_;
    my $maxNFeatures = 0;
    foreach my $cluster(@workingClusters){
        my @cluster = @$cluster;
        my $nFeatures = $cluster[$#cluster - 1];
        $maxNFeatures >= $nFeatures or $maxNFeatures = $nFeatures;
    }
    @$maxNClusters = ();
    @$remainingClusters = ();    
    foreach my $cluster(@workingClusters){  # split the remaining clusters for the next iteration
        my @cluster = @$cluster;
        my $nFeatures = $cluster[$#cluster - 1];
        my $array = $nFeatures == $maxNFeatures ? $maxNClusters : $remainingClusters;
        push @$array, $cluster;
    }
    my $nClusters = scalar(@$maxNClusters);
    my $nFeatures = $maxNFeatures * $nClusters;
    return ($maxNFeatures, $nClusters, $nFeatures);
}
sub runHotspotPermutation {  # run the collapsing permutation on the remaining features
    open my $permCollapsedH, ">", \my$permCollapsed or die $!;
    runPermutation($permCollapsedH, 1);
    close $permCollapsedH;
    loadClusters(\$permCollapsed, \my$nPermutedClusters, \my@permClusters);
    # relaxed = iterations with >=1 cluster with >= the current cluster feature count
    # strict  = iterations with >=N features in clusters with >= the current cluster feature count
    my (%relaxed, %strict);
    foreach my $cluster(@permClusters){
        my @cluster = @$cluster;
        my ($nFeatures, $iteration) = ($cluster[$#cluster - 1], $cluster[$#cluster]);
        $nFeatures >= $maxNFeatures and $relaxed{$iteration} += $nFeatures;
    }  
    foreach my $iteration(keys %relaxed){
         foreach my $nClusters (1..$nMaxClusters){
             $relaxed{$iteration} >= ($maxNFeatures * $nClusters) and $strict{$nClusters}{$iteration}++;
         }
    } 
    print join("\t", $nRemainingFeatures, $maxNFeatures, $nMaxClusters, $nMaxFeatures);
    foreach my $nClusters (1..$nMaxClusters){
        my $frequency = scalar(keys %{$strict{$nClusters}}) / $ENV{N_ITERATIONS};    
        print "\t$frequency";
    }
    print "\n";
}
sub excludeMaxNClusters {  # remove the highest ranking clusters from the previous iteration
    my $i = 0;
    my %remainingFeatures = map { $i++ => $_} @features;
    foreach my $cluster(@maxNClusters){
        push @{$exclusions{$$cluster[0]}{0}}, [$$cluster[1], $$cluster[2], 0, ['']];  # add cluster to exluded regions
        foreach my $i(keys %remainingFeatures){
            if ( $remainingFeatures{$i}[0] eq $$cluster[0] and
                 $remainingFeatures{$i}[1] >= $$cluster[1] and
                 $remainingFeatures{$i}[2] <= $$cluster[2] and
                 (!$ENV{STRANDED} or $remainingFeatures{$i}[5] eq $$cluster[5]) ) {
                delete $remainingFeatures{$i};  # delete the features within the cluster from the features list
            }
        }
    }
    @features = ();
    foreach my $i(keys %remainingFeatures){ push @features, $remainingFeatures{$i} }
    return scalar(@features); 
}
######################################################################

1;


