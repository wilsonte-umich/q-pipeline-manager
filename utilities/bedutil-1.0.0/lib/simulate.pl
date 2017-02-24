#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'simulate' takes feature BED lines on STDIN and streams them through
# the following bedutil functions:
#     permute | crossing | crosstab
# as a shortcut for assembling a complete simulation.  The final
# output is printed to STDOUT.  
#----------------------------------------------------------------------
# To view the options inherited from permute and crossing, use:
#     bedutil <command> --help
# Options must be passed as environment variables, not on the command line.
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 15, 9);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
);
%booleanOptions = map { $_ => 1} qw (
);
%defaultValues = (
);
@requiredOptions = qw (
);
my ($error, $feedback) = parseOptions($scriptName);
#######################################################################


#######################################################################
# assemble and execute the stream
#----------------------------------------------------------------------
my $permute =   getPerl("permute");
my $crossing =  getPerl("crossing");
my $crosstab =  getPerl("crosstab");
my $stream = "$permute | $crossing | $crosstab";
open my $outH, "|-", $stream or die "simulate.pl: could not open stream for writing: $!";
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

