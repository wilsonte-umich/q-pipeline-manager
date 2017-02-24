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
# handle output file backups
#-----------------------------------------------------------
sub backupInputs {
    my ($backupCommand) = @_;
    push @taskInputs, (
        $cgi->checkbox(-name=>'quiet',
                       -title=>"Do not show the names of files being (un)protected, backed up or restored"),           
    );
    $backupCommand eq 'restore' or return;
    push @taskInputs, (
        $cgi->checkbox(-name=>'force',
                       -title=>"Suppress warnings that duplicate jobs will be queued, files deleted, etc."), 
    );
}
sub getBackup {
    my ($backupCommand) = @_;
    my @options;
    setBooleanOptions(\@options, qw(quiet));
    $backupCommand eq 'restore' and setBooleanOptions(\@options, qw(force));
    my $qH = getQH([$backupCommand, @options]); 
    return parseQReturn($qH);  
}
############################################################

1;

