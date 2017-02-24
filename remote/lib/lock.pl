#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($params @taskInputs $updateMasterNames $debug);
our $isLocked = undef;
############################################################

############################################################
# q lock management
#-----------------------------------------------------------
sub setIsLocked {
    $$params{pageType} eq 'queue' or return;
    $$params{masterName} or return;
    defined $$params{isLocked} or return;
    my %isLocked = map { $_ => 1 } split(",", $$params{isLocked});
    $isLocked = $isLocked{$$params{masterName}};
}
sub toggleLock {
    $$params{pageType} eq 'queue' or return;
    $$params{toggleLock} and $updateMasterNames = 1;
}
sub lockInputs {
    my $state = $isLocked ? "locked" : "unlocked";
    push @taskInputs,
        "<p class=input>$$params{masterName} is currently $state</p>",
        "<input type=hidden name=toggleLock id=toggleLock>\n";
}
sub getLock {
    my $lockCommand = $isLocked ? "unlock" : "lock";
    my $qH = getQH([$lockCommand]); 
    return parseQReturn($qH); 
}
############################################################

1;

