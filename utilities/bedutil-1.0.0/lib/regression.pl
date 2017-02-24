#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'regression' takes feature BED lines on STDIN, streams them through 
# 'crossing' to collect a second score value, and then performs regression 
# analysis comparing the input scores (field 5) to the newly calculated scores.
# The correlation table is printed to STDOUT.  A plot of the simulation results
# is made and saved as a JPG image. 
#----------------------------------------------------------------------
# Options are:
#     JPG_FILE         name of the image file to create
#                      REQUIRED
#     X_LABEL          label for the X-axis, i.e. the name of the input values
#                      OPTIONAL [default: 'Input Value']
#     Y_LABEL          label for the Y-axis, i.e. the name of the calculated values
#                      OPTIONAL [default: 'Calculated Value']
#     MAX_INPUT_SCORE  maximum allowed value for the input scores
#                      any higher scores are set to $MAX_INPUT_SCORE prior to regression
#                      OPTIONAL [default: 0, no limit is enforced]
#     MAX_CALC_SCORE   label for the Y-axis, i.e. the name of the calculated values
#                      any higher scores are set to $MAX_CALC_SCORE prior to regression
#                      OPTIONAL [default: 0, no limit is enforced]
#----------------------------------------------------------------------
# To view the options inherited from crossing, use:
#     bedutil crossing --help
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil regression JPG_FILE=/path/to/my.jpg
# Options relevant to crossing must be passed as environment variables, 
# not on the command line.
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 35, 29);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    JPG_FILE  
    X_LABEL
    Y_LABEL
    MAX_INPUT_SCORE
    MAX_CALC_SCORE
);
%booleanOptions = map { $_ => 1} qw (
);
%defaultValues = (
    X_LABEL => 'Input Value',
    Y_LABEL => 'Calculated Value',
    MAX_INPUT_SCORE => 0,
    MAX_CALC_SCORE => 0
);
@requiredOptions = qw (
    JPG_FILE
);
my ($error, $feedback) = parseOptions($scriptName);
#######################################################################


#######################################################################
# assemble and execute the stream
#----------------------------------------------------------------------
my $crossing =  getPerl("crossing");
my $regression =  "Rscript $scriptDir/regression.R ".
  "$ENV{JPG_FILE} \"$ENV{X_LABEL}\" \"$ENV{Y_LABEL}\" $ENV{MAX_INPUT_SCORE} $ENV{MAX_CALC_SCORE}";
my $stream = "$crossing | $regression";
open my $outH, "|-", $stream or die "regression.pl: could not open stream for writing: $!";
while(my $line = <STDIN>){
    print $outH $line;
}
close $outH;
#----------------------------------------------------------------------
sub getPerl {
    my ($function) = @_;
    return "perl $scriptDir/$function.pl";
}
#######################################################################

1;

