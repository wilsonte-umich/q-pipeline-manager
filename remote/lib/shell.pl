#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @masterInputs @taskInputs @taskResults %project $scriptDir $isServer $debug);
############################################################

############################################################
# allow submission of commands directly to server shell
#-----------------------------------------------------------
sub shellInputs {
    @taskInputs = ();
    @masterInputs = ();    
    browseInputs();
    if($$params{browseDir} and !$$params{browseFileItem}){
        $$params{shellType} = 'shell';
        $$params{currentTask} = ""; 
        push @taskInputs, (
            $cgi->textarea(-name=>'shellCommand', 
                           -id=>'shellCommand', 
                           -class=>'shellCommand',),
            "<input type=hidden name=maskedShellCommand id=maskedShellCommand>",
            "<p class=input>Enter a command or multi-line script to be executed on server.  ".
            "Commands must be self-contained, e.g. variables are not preserved between calls.  ".
            "The starting directory is as specified above.</p>"
         )       
    } elsif($$params{browseFileItem}) {
        $$params{shellType} = 'file';
        $$params{currentTask} eq 'Edit' and getFileContents('browse');  
    }; 
    @taskInputs and unshift @taskInputs, '<hr>';
    push @taskInputs, (
        $cgi->hidden(-name=>'shellType', -id=>'shellType')
    )
}
#-----------------------------------------------------------
sub getShell {
    $$params{shellType} or return;
    if($$params{shellType} eq 'shell'){
        my $cmd = $$params{maskedShellCommand};
        $cmd or $cmd = "";
        $cmd =~ s/;\s*;/;/g;
        my $displayCmd = $cmd;
        $cmd = "cd $$params{browseDir} ; $cmd ";        
        unless($isServer){
            $cmd =~ s/\$/\\\$/g;
            $cmd =~ s/\"/\\\"/g;
            $cmd = '"'.$cmd.'"';
        }
        my $results = "EXECUTING:\n$displayCmd\n\nFROM:\n$$params{browseDir}\n\nRESULTS:\n".getSSH("$cmd 2>&1");
        open my $inH, "<", \$results;
        return parseQReturn($inH); 
    } elsif($$params{shellType} eq 'file'){
        getBrowse();
    }
}
############################################################

1;

