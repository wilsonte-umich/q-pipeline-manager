#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @taskInputs @taskResults $debug $scriptDir);
############################################################

############################################################
# q status web interface
#-----------------------------------------------------------
sub resubmitInputs {
    push @taskInputs, (
        $cgi->checkbox(-name=>'delete',
                       -title=>"Kill matching pending/running jobs when repeat job submissions are encountered"),
        $cgi->checkbox(-name=>'execute',
                       -title=>"Run target jobs immediately in shell instead of submitting with qsub"),
        $cgi->checkbox(-name=>'force',
                       -title=>"Suppress warnings that duplicate jobs will be queued, files deleted, etc."),                       
        $cgi->checkbox(-name=>'lock',
                       -title=>"Automatically lock masterFile after queuing jobs from it"),                       ,
        $cgi->checkbox(-name=>'no-chain',
                       -title=>"only apply command to --job; do not apply to jobs dependent on --job"),       
        "<div class=vSpacer></div>
        <table>",
        tablewrap(
            [$cgi->textfield(-name=>"job",
                             -id=>"job",
                             -title=>"Restrict command to specific jobID(s) (and sometimes its successors)",
	                         -class=>"wide"),
             "job"]),
        tablewrap(
            [$cgi->textfield(-name=>"depend",
                             -id=>"depend",
                             -title=>"Comma-delimited list of jobIDs upon which all new jobs should depend",
	                         -class=>"wide"),
             "depend"]),          
        "</table>",                           
    );
}
sub getResubmit {
    checkRequiredOptions(\@taskResults, qw(job));
    @taskResults and return;
    setBooleanOptions(\my@options, qw(delete execute force lock verbose no-chain));
    setValueOptions(\@options, qw(job depend));
    my $qH = getQH(['resubmit', @options]); 
    return parseQReturn($qH); 
}
############################################################

1;

