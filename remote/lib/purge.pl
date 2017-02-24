#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi @taskInputs);
############################################################

############################################################
# rollback and purge previously queued job files
#-----------------------------------------------------------
sub purgeInputs {
    push @taskInputs, (
        $cgi->checkbox(-name=>'force',
                       -title=>"Suppress warnings that duplicate jobs will be queued, files deleted, etc."), 
    );
}
sub getPurge {
    my ($purgeCommand) = @_;
    setBooleanOptions(\my@options, qw(force));
    my $qH = getQH([$purgeCommand, @options]); 
    return parseQReturn($qH);  
}
############################################################

1;

