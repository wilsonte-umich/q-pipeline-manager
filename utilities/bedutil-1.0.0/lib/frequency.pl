#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'frequency' takes feature BED lines on STDIN as input that have had an
# evaluation score added as the last column, and a simulation iteration
# added as the penultimate column.  For each iteration, an aggregate feature
# score is calculated.  A frequency distribution of these aggregrate scores
# across all iterations is constructed and printed to STDOUT, along with
# the offset control (input iteration -2), if present, in the first row,
# the  actual aggregrate value in the second row (input iteration -1), and a
# header row.  If instructed, a plot image of the simulation is also made.
#----------------------------------------------------------------------
# Options are:
#     AGGREGATE_TYPE  how to aggregate the feature scores within each iteration
#                     allowed values are enumerated below
#                     OPTIONAL [default: AVERAGE]
#     N_ITERATIONS    the number of permutation iterations represented by the data
#                     OPTIONAL [default: determined from input data]
#     JPG_FILE        file path of jpg graph to be created
#                     OPTIONAL [default: no graph is created]
#     X_LABEL         label for the X-axis, i.e. the name of the evaluation score
#                     OPTIONAL [default: 'Aggregate Score']
#     SUPPRESS_NULL   boolean instructing to ignore input features with score of 'NULL'
#                     if TRUE, aggregate values will be determined from only non-NULL scores
#                     OPTIONAL [default: FALSE, 'NULL' scores are interpreted as 0]
#----------------------------------------------------------------------
# Aggregate types are:
#     AVERAGE         the average of all scores in an iteration
#     SUM             the sum of all scores in an iteration
#     COUNT           the number of boolean-positive iteration scores
#     MEDIAN          the median of all scores in an iteration
#     MIN             the minimum iteration score value
#     MAX             the maximum iteration score value
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>
# with no spaces, for example:
#     bedutil frequency X_LABEL="My Data"
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 41, 35);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    AGGREGATE_TYPE
    N_ITERATIONS
    JPG_FILE
    X_LABEL
    SUPPRESS_NULL
);
%booleanOptions = map { $_ => 1} qw (
    SUPPRESS_NULL
);
%defaultValues = (
    AGGREGATE_TYPE => 'AVERAGE',
    X_LABEL => 'Aggregate Score'
);
@requiredOptions = qw (
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
$ENV{IS_CHILD} and return 1;
#######################################################################


#######################################################################
# collect BED features from STDIN
#----------------------------------------------------------------------
loadBedFile(undef, $error, \my$nFeatures, \my@features);
$nFeatures or exit;
print STDERR "$feedback: $nFeatures features provided on STDIN\n";
my $out = "cat";
$ENV{JPG_FILE} and $out = "Rscript $scriptDir/frequency.R $ENV{JPG_FILE} \"$ENV{X_LABEL}\"";
open my $outH, "|-", $out or die "frequency.pl: could not open stream for writing: $!";
my (%counts, %sums, %positives, %scores, %fixedAggs, %iterAggs, %iterHist);
my ($nLow, $nHigh) = (0, 0);
printFrequencies($outH, \@features);
close $outH;
#######################################################################


#######################################################################
# worker subs
#----------------------------------------------------------------------
sub printFrequencies {
    my ($outH, $features) = @_;
    collectIterationScores($features);
    getAggregateScores();
    my $nIterations = getNIterations();
    print $outH join("\t", qw(aggregateScore frequency fractionAtOrBelow fractionAtOrAbove)), "\n";
    print $outH join("\t", $fixedAggs{-2}, -2, -2, -2), "\n";
    print $outH join("\t", $fixedAggs{-1}, -1, $nLow / $nIterations, $nHigh / $nIterations), "\n";
    foreach my $roundedAgg(sort {$a <=> $b} keys %iterHist){
        my $frequency = $iterHist{$roundedAgg} / $nIterations;
        my ($nLow, $nHigh) = (0, 0);
        foreach my $iterAgg(keys %iterAggs){
            $iterAgg <= $roundedAgg and $nLow += $iterAggs{$iterAgg};
            $iterAgg >= $roundedAgg and $nHigh += $iterAggs{$iterAgg};
        }
        print $outH join("\t", $roundedAgg, $frequency, $nLow / $nIterations, $nHigh / $nIterations), "\n";
    }
}
sub collectIterationScores {  # make initial aggregating calculations on the input iteration sets
    my ($features) = @_;
    foreach my $feature(@$features){
        my $score = pop @$feature;
        if($score eq 'NULL'){
            $ENV{SUPPRESS_NULL} and next;
            $score = 0;
        }
        my $iteration = pop @$feature;
        $counts{$iteration}++;
        $sums{$iteration} += $score;
        $positives{$iteration} += $score ? 1 : 0;  # interpret score as a boolean
        push @{$scores{$iteration}}, $score;
    }
}
sub getAggregateScores {  # make the final aggegrate calculations per iteration
    my ($minIterAgg, $maxIterAgg, @iterAggs);
    foreach my $iteration(sort {$a <=> $b} keys %counts){
        my $aggregate = getAggregateScore($iteration);
        if($iteration < 0){
            $fixedAggs{$iteration} = $aggregate;
            my $iterName = $iteration == -1 ? "actual" : "offset";
            print STDERR "$feedback: $aggregate\t$iterName aggregated score\n";
        } else {
            defined $fixedAggs{-1} or die "$error: actual data must be present with iteration == -1\n";
            (defined $minIterAgg and $minIterAgg <= $aggregate) or $minIterAgg = $aggregate;
            (defined $maxIterAgg and $maxIterAgg >= $aggregate) or $maxIterAgg = $aggregate;
            $iterAggs{$aggregate}++;
            $aggregate <= $fixedAggs{-1} and $nLow++;
            $aggregate >= $fixedAggs{-1} and $nHigh++;
        }
    }
    print STDERR "$feedback: $minIterAgg\tlowest aggregated iteration score\n";
    print STDERR "$feedback: $maxIterAgg\thighest aggregated iteration score\n";
    my $binSize = ($maxIterAgg - $minIterAgg)/50;
    $binSize or die "$error: no range of iteration aggregate values was observed for binning\n";
    foreach my $iterAgg(keys %iterAggs){
        my $roundedAgg = int($iterAgg / $binSize + 0.5) * $binSize;
        $iterHist{$roundedAgg} += $iterAggs{$iterAgg};
    }
}
sub getAggregateScore {  # return the requested iteration aggregate score
    my ($iteration) = @_;
    if($ENV{AGGREGATE_TYPE} eq 'AVERAGE'){
        return $sums{$iteration} / $counts{$iteration};
    } elsif($ENV{AGGREGATE_TYPE} eq 'SUM'){
        return $sums{$iteration};
    } elsif($ENV{AGGREGATE_TYPE} eq 'COUNT'){
        return $positives{$iteration};
    } elsif($ENV{AGGREGATE_TYPE} eq 'MEDIAN'){
        my @scores = sort {$a <=> $b} @{$scores{$iteration}};
        my $i = int($counts{$iteration}/2+0.5);
        return $scores[$i];
    } elsif($ENV{AGGREGATE_TYPE} eq 'MIN'){
        my @scores = sort {$a <=> $b} @{$scores{$iteration}};
        return $scores[0];
    } elsif($ENV{AGGREGATE_TYPE} eq 'MAX'){
        my @scores = sort {$b <=> $a} @{$scores{$iteration}};
        return $scores[0];
    } else {
        die "$error: unknown aggregate type: $ENV{AGGREGATE_TYPE}\n";
    }
}
sub getNIterations{  # collect the number of iterations
    delete $counts{-2};
    delete $counts{-1};
    my $nIterations = $ENV{N_ITERATIONS} ? $ENV{N_ITERATIONS} : scalar(keys %counts);
    $nIterations or die "$error: N_ITERATIONS was 0 or absent\n";
    return $nIterations;
}
#######################################################################

1;

