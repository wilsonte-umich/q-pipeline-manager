#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($libDir $cgi $params %project $isSSH 
            $isCacheRequest $isServer $debug);
our (@taskButtons, @subTaskButtons, @taskInputs, @taskResults);
############################################################

############################################################
# determine desired task, inputs and results
#----------------------------------------------------------
sub taskButtons {
    @taskButtons = ();
    if($$params{pageType} eq 'queue'){
        queueTaskButtons();
    } elsif($$params{pageType} eq 'shell'){
        browseTaskButtons();
    }
    @taskButtons and unshift @taskButtons, '<hr>';    
    push @taskButtons, $cgi->hidden(-name=>"currentTask", -id=>"currentTask"),
}
sub taskButton {
    my ($label, $task, $title) = @_;
    my $id = $task."Button";
    "<input type=button class=activeButton id=$id value=$label ".
            "onClick=setCurrentTask('$task') title=\"$title\">&nbsp";
}
#-----------------------------------------------------------
sub subTaskButtons {
    @subTaskButtons = ();
    $isSSH or return;
    if($$params{pageType} eq 'queue'){
        queueSubTaskButtons();
    } elsif($$params{pageType} eq 'browse'){
        $$params{currentSubTask} = "";
    } elsif($$params{pageType} eq 'shell'){
        $$params{currentSubTask} = "";
    }
    @subTaskButtons and unshift @subTaskButtons, '<hr>';    
    push @subTaskButtons, $cgi->hidden(-name=>"currentSubTask", -id=>"currentSubTask");
}
sub subTaskButton {
    my ($label, $subTask, $title) = @_;
    my $id = $subTask."Button";
    "<input type=button class=activeButton id=$id value=$label ".
            "onClick=setCurrentSubTask('$subTask') title=\"$title\">&nbsp";
}
sub subTaskSilencedButton {
    my ($label, $subTask, $title) = @_;
    my $id = $subTask."Button";
    "<input type=button class=silencedButton id=$id value=$label ".
            "title=\"$title\">&nbsp";
}
#-----------------------------------------------------------
sub taskInputs {
    @taskInputs = ();
    if($$params{pageType} eq 'queue'){
        queueInputs();
    } elsif($$params{pageType} eq 'shell'){
        shellInputs();
    }     
}
#-----------------------------------------------------------
sub taskResults {
    $$params{editingFile} and return;
    $isServer and @taskResults = ();
    $$params{executing} or return;
    $isCacheRequest and return;  # using cached results
    @taskResults = ();
    if($$params{pageType} eq 'queue'){
        getQueue();
    } elsif($$params{pageType} eq 'shell'){
        getShell();
    }  
}
############################################################

1;

