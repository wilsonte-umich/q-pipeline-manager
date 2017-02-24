#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'relate' takes feature BED lines that have had a dependent evaluation score 
# added as the last column, a simulation iteration added as the penultimate
# column, and that bear an independent feature score in BED column 5.
# A plot is constructed correlating the dependent score to the independent
# score for the actual data points, showing ranges for the simulation. 
# A correlation coefficient is calculated for the actual data.  Data are
# provided on STDIN.  The output is an image file.
#----------------------------------------------------------------------
# Options are:
#     JPG_FILE       file path of jpg graph to be created
#                    REQUIRED
#     INDEPENDENT    name/label of the independent variable, plotted on the X axis
#                    OPTIONAL [default: 'Independent']
#     DEPENDENT      name/label of the dependent variable, plotted on the Y axis
#                    OPTIONAL [default: 'Dependent']
#     SUPPRESS_SIM   boolean indicating whether to suppress the simulation ranges
#                    OPTIONAL [default: FALSE, the simulation ranges are plotted]
#     X_LOG          boolean indicating whether to log10 transform the independent scores
#                    OPTIONAL [default: FALSE]
#     X_MIN          lower limit of the X axis
#                    data below this limit are plotted and aggregated at the lower limit
#                    OPTIONAL [default: 0, or smallest value if X_LOG=TRUE]
#     X_MAX          upper limit of the X axis
#                    data above this limit are plotted and aggregated at the upper limit
#                    OPTIONAL [default: largest value in data]
#     Y_LOG          boolean indicating whether to log10 transform the dependent scores
#                    OPTIONAL [default: FALSE]
#     Y_MIN          lower limit of the Y axis
#                    data below this limit are plotted at the lower limit
#                    OPTIONAL [default: 0, or smallest value if Y_LOG=TRUE]
#     Y_MAX          upper limit of the Y axis
#                    data above this limit are plotted at the upper limit
#                    OPTIONAL [default: largest value in data]
#     ROUND_DIGITS   number of decimal places to round dependent values by prior to aggregation
#                    rounding is applied after log is applied, if applicable
#                    OPTIONAL [default: no rounding is applied, dependent values are used as is]
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil relate JPG_FILE=my.jpg INDEPENDENT="My Value"
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 47, 41);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    JPG_FILE
    INDEPENDENT    
    DEPENDENT
    SUPPRESS_SIM
    X_LOG
    X_MIN
    X_MAX
    Y_LOG
    Y_MIN
    Y_MAX
    ROUND_DIGITS
);
%booleanOptions = map { $_ => 1} qw (
    SUPPRESS_SIM
    X_LOG
    Y_LOG
);
%defaultValues = (
    INDEPENDENT => 'Independent',
    DEPENDENT => 'Dependent',
    X_MIN => '',
    X_MAX => '',
    Y_MIN => '',
    Y_MAX => '',
    ROUND_DIGITS => '',
);
@requiredOptions = qw (
    JPG_FILE
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
$ENV{IS_CHILD} and return 1;
#######################################################################


#######################################################################
# collect BED features from STDIN
#----------------------------------------------------------------------
print "relating \"$ENV{INDEPENDENT}\" to \"$ENV{DEPENDENT}\"\n\n";
my $out = "Rscript $scriptDir/relate.R";   
open my $outH, "|-", $out or die "relate.pl: could not open stream for writing: $!";
while(my $line = <STDIN>){
    print $outH $line;
}
close $outH;
print "\ndone\n";
#######################################################################

1;

