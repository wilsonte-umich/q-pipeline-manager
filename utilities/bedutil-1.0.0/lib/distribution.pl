#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'distribution' takes feature BED lines on STDIN as input that have had an
# evaluation score added as the last column, and a simulation iteration
# added as the penultimate column. Four frequency distributions are plotted:
#     unpaired_all  set 1 = actual scores, set 2 = all individual iteration scores
#     unpaired_agg  set 1 = actual scores, set 2 = iteration scores aggregated by feature
#     paired_all    difference between actual scores and all paired individual iteration scores
#     paired_agg    difference between actual scores and iteration scores aggregated by feature
# The BED name field, field 4, is used as the feature unique identifier
# for pairing actual to simulated features.
#----------------------------------------------------------------------
# Options are:
#     IMAGE_PREFIX     path prefix of the various image files that will be created
#                      REQUIRED
#     X_LABEL          label for the X-axis, i.e. the name of the evaluation score
#                      OPTIONAL [default: 'Score']
#     SUPPRESS_NULL    boolean instructing to ignore input features with score of 'NULL'
#                      OPTIONAL [default: FALSE, 'NULL' scores are interpreted as 0]
#     N_BINS           number of bins in the histogram
#                      OPTIONAL [default: 50]
#     AGGREGATE_TYPE   how to aggregate the iteration scores for each input feature
#                      allowed values are enumerated below
#                      OPTIONAL [default: MEDIAN]
#----------------------------------------------------------------------
# Aggregate types are:
#     AVERAGE          the average of all iteration scores for an input feature
#     MEDIAN           the median of all iteration scores for an input feature
#     MIN              the minimum iteration score for the feature
#     MAX              the maximum iteration score for the feature
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>
# with no spaces, for example:
#     bedutil distribution X_LABEL="My Score" SUPPRESS_NULL
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 38, 32);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    IMAGE_PREFIX
    X_LABEL
    SUPPRESS_NULL
    N_BINS
    AGGREGATE_TYPE
);
%booleanOptions = map { $_ => 1} qw (
    SUPPRESS_NULL
);
%defaultValues = (
    X_LABEL => 'Score',
    N_BINS => 50,
    AGGREGATE_TYPE => 'MEDIAN'
);
@requiredOptions = qw (
    IMAGE_PREFIX
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
$ENV{IS_CHILD} and return 1;
my %aggFunctions = (
    AVERAGE => 'mean',
    MEDIAN  => 'median',
    MIN     => 'min',
    MAX     => 'max',
);
$ENV{AGGREGATE_FUNCTION} = $aggFunctions{$ENV{AGGREGATE_TYPE}};
#######################################################################


#######################################################################
# assemble and execute the stream
#----------------------------------------------------------------------
my $stream = "perl $scriptDir/crosstab.pl | Rscript $scriptDir/distribution.R";  # crosstab the permute+crossing output
open my $outH, "|-", $stream or die "distribution.pl: could not open stream for writing: $!";
while(my $line = <STDIN>){
    print $outH $line;
}
close $outH;
#######################################################################

1;

