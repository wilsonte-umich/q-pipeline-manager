#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'cdf' takes feature BED lines on STDIN containg both actual and iteration
# feature scores as formatted by the crosstab command. The estimated 
# cumulative distribution function (CDF) of the iterations of the smallest, 
# median, and largest input features are plotted to visualize the nature of 
# the score distribution and how it changes with feature size. Reference  
# lines are plotted corresponding to the normal CDF at the mean and standard 
# deviation of each of the three iteration sets. A plot is created for 
# unmodified iteration values as well as log10 transformed values. 
# A separate plot relates score rank to quantiles. Unmodified input lines are 
# repeated to STDOUT for analysis streaming.
#----------------------------------------------------------------------
# Options are:
#     IMAGE_PREFIX     path prefix of the image files that will be created
#                      suffixed with ".unmodified_cdf.jpg/pdf", ".log10_cdf.jpg/pdf" and "rank_quantile.jpg/pdf"
#                      REQUIRED
#     SIM_LABEL        label for the X-axis, i.e. the name of the simulation score
#                      OPTIONAL [default: 'Score']
#     SUPPRESS_NULL    boolean instructing to ignore input features with score of 'NULL'
#                      OPTIONAL [default: FALSE, 'NULL' scores are interpreted as 0]           
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>
# with no spaces, for example:
#     bedutil cdf SUPPRESS_NULL > /dev/null
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 30, 24);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    IMAGE_PREFIX
    SIM_LABEL
    SUPPRESS_NULL
);
%booleanOptions = map { $_ => 1} qw (
    SUPPRESS_NULL
);
%defaultValues = (
    SIM_LABEL => 'Score',
);
@requiredOptions = qw (
    IMAGE_PREFIX
);
my ($error, $feedback) = parseOptions($scriptName, $scriptDir);
#######################################################################


#######################################################################
# assemble and execute the stream
#----------------------------------------------------------------------
print STDERR "bedutil cdf\n";
my $stream = "Rscript $scriptDir/cdf.R";  # crosstab the permute+crossing output
open my $outH, "|-", $stream or die "cdf.pl: could not open stream for writing: $!";
while(my $line = <STDIN>){
    print $outH $line;
}
close $outH;
#######################################################################

1;

