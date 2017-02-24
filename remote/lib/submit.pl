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
sub submitInputs {
    push @taskInputs, (
        $cgi->checkbox(-name=>'delete',
                       -title=>"Kill matching pending/running jobs when repeat job submissions are encountered"),
        $cgi->checkbox(-name=>'execute',
                       -title=>"Run target jobs immediately in shell instead of submitting with qsub"),
        $cgi->checkbox(-name=>'force',
                       -title=>"Suppress warnings that duplicate jobs will be queued, files deleted, etc."),                       
        $cgi->checkbox(-name=>'lock',
                       -title=>"Automatically lock masterFile after queuing jobs from it"),                       
        $cgi->checkbox(-name=>'verbose',
                       -title=>"Report all commands acted on from masterFile (extremely verbose)"),
        "<div class=vSpacer></div>
        <table>",
        tablewrap(
            [$cgi->textfield(-name=>"depend",
                             -id=>"depend",
                             -title=>"Comma-delimited list of jobIDs upon which all new jobs should depend",
	                         -class=>"wide"),
             "depend"],
        ),           
        "</table>",                           
    );
}
sub getSubmit {
    my ($submitCommand) = @_;
    setBooleanOptions(\my@options, qw(delete execute force lock verbose));
    setValueOptions(\@options, qw(depend));
    my $qH = getQH([$submitCommand, @options]); 
    return parseQReturn($qH); 
}
############################################################

1;

