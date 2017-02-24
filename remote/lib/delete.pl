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
# q delete web interface
#-----------------------------------------------------------
sub deleteInputs {
    push @taskInputs, (
        $cgi->checkbox(-name=>'force',
                       -title=>"Suppress warnings that duplicate jobs will be queued, files deleted, etc."),                       
        "<div class=vSpacer></div>
        <table>",
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
sub getDelete {
    checkRequiredOptions(\@taskResults, qw(job));
    @taskResults and return;
    setBooleanOptions(\my@options, qw(force));
    setValueOptions(\@options, qw(job));
    my $qH = getQH(['delete', @options]); 
    return parseQReturn($qH);  
}
############################################################

1;

