#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @taskInputs @taskResults);
############################################################

############################################################
# q clear web interface
#-----------------------------------------------------------
sub clearInputs {
    push @taskInputs, (
        "<table>",
        tablewrap(
            [$cgi->textfield(-name=>"job",
                             -id=>"job",
                             -title=>"Restrict command to specific jobID(s) (and sometimes its successors)",
	                         -class=>"wide"),
             "job"],
        ),           
        "</table>",
    );
}
sub getClear {
    checkRequiredOptions(\@taskResults, qw(job));
    @taskResults and return;
    setValueOptions(\my@options, qw(job));
    my $qH = getQH(['clear', @options]); 
    return parseQReturn($qH);  
}
############################################################

1;

