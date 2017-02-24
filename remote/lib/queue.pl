#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params $isSSH @taskButtons @subTaskButtons 
            @taskInputs $isLocked);
my %subTasks = (_edit =>    {master=>1,template=>1,help=>1,shell=>1},
                _submit =>  {submit=>1,extend=>1,lock=>1,rollback=>1,purge=>1},
                _monitor => {status=>1,clear=>1,delete=>1,resubmit=>1},
                _protect => {protect=>1,unprotect=>1,backup=>1,restore=>1},
                _publish => {publish=>1} );
my %defaultSubTasks = (_edit=>'master',_submit=>'submit',_monitor=>'status',_protect=>'protect',_publish=>'publish');
############################################################

############################################################
# switchboard to inputs and results for q itself
#-----------------------------------------------------------
sub queueTaskButtons {
    unless($$params{masterName}){ 
        $$params{currentTask} = ""; 
        $$params{currentSubTask} = ""; 
        $$params{isStatus} = ""; 
        return;
    }
    $isSSH or return;
    $$params{masterNameChanged} and delete $$params{currentTask};
    $$params{currentTask} or $$params{currentTask} = '_edit';
    push @taskButtons, (
        taskButton("Edit",  
                   "_edit",
                   "View/edit the master, README or template files"),
        taskButton("Submit",  
                   "_submit",
                   "Check the syntax of a master file and submit its jobs to the scheduler"),
        taskButton("Monitor", 
                   "_monitor",
                   "Monitor or modify jobs previously queued by a master file"),
        taskButton("Protect", 
                   "_protect",
                   "Write protect and back up pipeline output files"),
        taskButton("Publish", 
                   "_publish",
                   "Create and view publication reports for a master file"),
        $cgi->hidden(-name=>"isStatus", -id=>"isStatus"),
    );
}
#-----------------------------------------------------------
sub queueSubTaskButtons {
    my $ct = $$params{currentTask}; 
    $ct or return;  
    $$params{masterNameChanged} and delete $$params{currentSubTask};
    $$params{currentSubTask} or $$params{currentSubTask} = $defaultSubTasks{$ct};
    $subTasks{$ct}{$$params{currentSubTask}} or $$params{currentSubTask} = $defaultSubTasks{$ct};
    my $cst = $$params{currentSubTask};       
    if($ct eq '_edit'){
        push @subTaskButtons, (
            subTaskButton("Master",  
                          "master",
                          "View/edit the selected master file"),
            subTaskButton("README",  
                          "help",
                          "View/edit the master class README file"),
            subTaskButton("Template",  
                          "template",
                          "View/edit the master class template"),
        );
    } elsif($ct eq '_submit') {
        my $lockLabel = $isLocked ? "Unlock" : "Lock";
        my $lockTitle = $isLocked ? 
                        "Remove the protective marker set by Lock" : 
                        "Set a protective marker that prevents job actions from masterFile";
        push @subTaskButtons, (
            subTaskButton("Submit",  
                          "submit",
                          "Submit all jobs to the job scheduler"),
            subTaskButton("Extend",  
                          "extend",
                          "Submit previously unsatisfied jobs to the job scheduler"),
            subTaskButton($lockLabel,  
                          "lock",
                          $lockTitle),
            subTaskButton("Rollback",  
                          "rollback",
                          "Revert pipeline to the most recently archived status file"),
            subTaskButton("Purge",  
                          "purge",
                          "Remove all status, script and log files created by masterFile")                 
        );
    } elsif($ct eq '_monitor') {
        push @subTaskButtons, (
            subTaskButton("Status",  
                          "status",
                          "Update and show the status of queued jobs"),
            subTaskButton("Clear",  
                     "clear",
                     "Clear error states from specified jobs"),
            subTaskButton("Delete",  
                     "delete",
                     "Kill (qdel) incomplete jobs"),
            subTaskButton("Resubmit",  
                     "resubmit",
                     "Submit identical copies of previously queued jobs"),
        ); 
    } elsif($ct eq '_protect') {
        push @subTaskButtons, (
            subTaskButton("Protect",  
                          "protect",
                          "Write-protect files identified as 'protect <fileGlob>' in masterFile"),
            subTaskButton("Unprotect",  
                          "unprotect",
                          "Remove file write-protection for --who (chmod <--who>+w)"),
            subTaskButton("Backup",  
                          "backup",
                          "Create a copy of directories identified as 'backup <directory>' in masterFile"),
            subTaskButton("Restore",  
                          "restore",
                          "Restore from backup directories identified as 'backup <directory>' in masterFile"), 
        ); 
    } elsif($ct eq '_publish') {
        push @subTaskButtons, (
            subTaskButton("Publish",  
                          "publish",
                          "Create a distribution of instruction, script, status, and log files")
        ); 
    }
}
#-----------------------------------------------------------
sub queueInputs {
    my $ct = $$params{currentTask};
    my $cst = $$params{currentSubTask}; 
    ($ct and $cst) or return;
    if($ct eq '_edit'){
        if($cst eq 'master') {
            editInputs($cst);
        } elsif($cst eq 'template'){
            editInputs($cst);
        } elsif($cst eq 'help'){
            editInputs($cst);  
        }
    } elsif($ct eq '_submit'){
        if($cst eq 'submit') {
            submitInputs();
        } elsif($cst eq 'extend'){
            submitInputs();
        } elsif($cst eq 'lock'){
            lockInputs();
        } elsif($cst eq 'rollback'){
            purgeInputs();
        } elsif($cst eq 'purge'){
            purgeInputs();
        }  
    } elsif($ct eq '_monitor'){
        if($cst eq 'status') {
            statusInputs();
        } elsif($cst eq 'clear'){
            clearInputs();
        } elsif($cst eq 'delete'){
            deleteInputs();
        } elsif($cst eq 'resubmit'){
            resubmitInputs();
        }
    } elsif($ct eq '_protect'){
        if($cst eq 'protect') {
            protectInputs($cst);
        } elsif($cst eq 'unprotect'){
            protectInputs($cst);
        } elsif($cst eq 'backup'){
            backupInputs($cst);
        } elsif($cst eq 'restore'){
            backupInputs($cst);
        }
    } elsif($ct eq '_publish'){
        if($cst eq 'publish') {
            publishInputs($cst);
        } 
    }
    @taskInputs and unshift @taskInputs, "<hr><input type=hidden name=getReport id=getReport >";
}
#---------------------------------------------------------
sub getQueue { # q task results
    my $ct = $$params{currentTask}; 
    my $cst = $$params{currentSubTask}; 
    if($$params{getReport}){  # redirect jobID click to report page
        getReport();      
    } elsif($ct eq '_edit'){
        if($cst eq 'master') {
            saveFile($cst);
        } elsif($cst eq 'template'){
            saveFile($cst);
        } elsif($cst eq 'help'){
            saveFile($cst);  
        }
    } elsif($ct eq '_submit'){
        if($cst eq 'submit') {
            getSubmit($cst);
        } elsif($cst eq 'extend'){
            getSubmit($cst);
        } elsif($cst eq 'lock'){
            getLock();
        } elsif($cst eq 'rollback'){
            getPurge($cst);
        } elsif($cst eq 'purge'){
            getPurge($cst);
        }      
    } elsif($ct eq '_monitor'){
        if($cst eq 'status') {
            getStatus();
        } elsif($cst eq 'clear'){
            getClear();
        } elsif($cst eq 'delete'){
            getDelete();
        } elsif($cst eq 'resubmit'){
            getResubmit();
        } 
    } elsif($ct eq '_protect'){
        if($cst eq 'protect') {
            getProtect($cst);
        } elsif($cst eq 'unprotect'){
            getProtect($cst);
        } elsif($cst eq 'backup'){
            getBackup($cst);
        } elsif($cst eq 'restore'){
            getBackup($cst);
        } 
    } elsif($ct eq '_publish'){
        if($cst eq 'publish') {
            getPublish($cst);
        } 
    }
}
############################################################

1;

