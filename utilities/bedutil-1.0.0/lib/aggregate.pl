#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'aggregate' takes feature BED lines on STDIN containg both actual and iteration
# feature scores as formatted by the crosstab command and estimates enrichment
# p-values by two methods. 1) The distribution of aggregated feature scores for
# each iteration is plotted and compared to the actual aggregate feature score,
# printing to STDERR the fraction of iterations higher or lower than the actual
# aggregated score and a one-sided p-value estimated from a Gaussian curve fit
# 2) The distribution of the quantiles of the actual scores for each feature
# relative to its iterations is plotted and subjected to one-sample two-sided
# Sign and Wilcoxon tests. p-values are only calculated if the underlying
# score distributions permit: the Gaussian estimate requires that the iteration
# aggregate scores are normally distributed; the sign test requires that a true
# median iteration score can be defined for all features; the Wilcoxon test
# requires that the iteration scores are continuously distributed for all
# features. Unmodified input lines are repeated to STDOUT for analysis streaming.
#----------------------------------------------------------------------
# Options are:
#     IMAGE_PREFIX     path prefix of the image files that will be created
#                      suffixed with ".by_iteration.jpg/pdf" and ".by_feature.jpg/pdf"
#                      REQUIRED
#     SIM_LABEL        label for the simulation score type
#                      OPTIONAL [default: 'Simulation score']
#     SUPPRESS_NULL    boolean instructing to ignore input features with simulation score of 'NULL'
#                      OPTIONAL [default: FALSE, 'NULL' scores are interpreted as 0]
#     SIM_LOG          boolean indicating to take the log10 of simulation scores before analysis
#                      OPTIONAL [default: FALSE]
#     AGGREGATE_TYPE   how to aggregate the feature simulation scores within each iteration
#                      allowed values are enumerated below
#                      OPTIONAL [default: MEDIAN]
#     N_BINS           number of bins in the simulation score histograms
#                      OPTIONAL [default: 50]
#----------------------------------------------------------------------
# Aggregate types are:
#     AVERAGE         the average of all simulation scores
#     WEIGHTED        the average of all simulation scores weighted by feature length
#     SUM             the sum of all simulation scores
#     COUNT           the number of boolean-positive simulation scores
#     MEDIAN          the median of all scores
#     MIN             the minimum simulation score value
#     MAX             the maximum simulation score value
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>
# with no spaces, for example:
#     bedutil aggregate AGGREGATE_TYPE=AVERAGE SUPPRESS_NULL > /dev/null
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 50, 44);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    IMAGE_PREFIX
    SIM_LABEL
    SUPPRESS_NULL
    SIM_LOG
    AGGREGATE_TYPE
    N_BINS
);
%booleanOptions = map { $_ => 1} qw (
    SUPPRESS_NULL
    SIM_LOG
);
%defaultValues = (
    SIM_LABEL => "Simulation score",
    AGGREGATE_TYPE => "MEDIAN",
    N_BINS => 50
);
@requiredOptions = qw (
    IMAGE_PREFIX
);
my ($error, $feedback) = parseOptions($scriptName, $scriptDir);
#######################################################################


#######################################################################
# assemble and execute the stream
#----------------------------------------------------------------------
print STDERR "bedutil aggregate\n";
my $stream = "Rscript $scriptDir/aggregate.R";
open my $outH, "|-", $stream or die "aggregate.pl: could not open stream for writing: $!";
while(my $line = <STDIN>){
    print $outH $line;
}
close $outH;
#######################################################################

1;

