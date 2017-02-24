#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'normtest' takes feature BED lines on STDIN containg both actual and iteration
# feature scores as formatted by the crosstab command. The Shapiro-Wilk
# normality test is applied to the iteration values of each feature. A plot
# is constructed showing the normality test p-value as a function of feature 
# length, for unmodified iteration values as well as log10 transformed values. 
# Unmodified input lines are repeated to STDOUT for analysis streaming.
#----------------------------------------------------------------------
# Options are:
#     IMAGE_PREFIX     path prefix of the image file that will be created
#                      suffixed with ".normtest.jpg/pdf"
#                      REQUIRED
#     SIM_LABEL        title for the plot, e.g. the name of the simulation score
#                      OPTIONAL [default: '']
#     SUPPRESS_NULL    boolean instructing to ignore input features with score of 'NULL'
#                      OPTIONAL [default: FALSE, 'NULL' scores are interpreted as 0]            
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>
# with no spaces, for example:
#     bedutil normtest SUPPRESS_NULL > /dev/null
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 26, 20);
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
    SIM_LABEL => ""
);
@requiredOptions = qw (
    IMAGE_PREFIX
);
my ($error, $feedback) = parseOptions($scriptName, $scriptDir);
#######################################################################


#######################################################################
# assemble and execute the stream
#----------------------------------------------------------------------
print STDERR "bedutil normtest\n";
my $stream = "Rscript $scriptDir/normtest.R";  
open my $outH, "|-", $stream or die "normtest.pl: could not open stream for writing: $!";
while(my $line = <STDIN>){
    print $outH $line;
}
close $outH;
#######################################################################

1;

