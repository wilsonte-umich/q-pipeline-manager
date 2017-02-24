#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi @taskInputs @taskResults);
############################################################

############################################################
# handle output file write protection
#-----------------------------------------------------------
sub protectInputs {
    my ($protectCommand) = @_;
    push @taskInputs, (
        $cgi->checkbox(-name=>'quiet',
                       -title=>"Do not show the names of files being (un)protected, backed up or restored"),           
    );
    $protectCommand eq 'unprotect' or return;
    push @taskInputs, (
        "<table>",
        tablewrap(
            [$cgi->textfield(-name=>"who",
                             -id=>"who",
                             -title=>"List of classes to unprotect, consistent with chmod <--who>+w (e.g. a, u, or g)",
	                         -class=>"wide",
	                         -value=>"ug"),
             "who"],
        ),           
        "</table>",
    );
    
}
sub getProtect {
    my ($protectCommand) = @_;
    my @options;
    if($protectCommand eq 'unprotect'){
        checkRequiredOptions(\@taskResults, qw(who));
        @taskResults and return;
        setValueOptions(\@options, qw(who));
    }
    setBooleanOptions(\@options, qw(quiet));
    my $qH = getQH([$protectCommand, @options]); 
    return parseQReturn($qH);  
}
############################################################

1;

