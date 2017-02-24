#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($imagesDir %config $cgi $params %project @taskInputs 
            @masterInputs @taskButtons @taskResults $isSSH $isServer $debug);
my @fileTasks = qw(head tail cat stat Edit);
my %fileTasks = (
    head => 'Show the first lines of a file',
    tail => 'Show the last lines of a file',
    cat => 'Show the entire contents of a file',
    stat => 'Show basic statistics for file',
    Edit => 'Edit a file',
);
my %imageExtensions = map { $_ => $_ } qw(jpg png gif bmp tif);
############################################################

############################################################
# file system browsing
#-----------------------------------------------------------
sub browseInputs {
    my $startDir = $project{directory};
    $startDir or $startDir = "";
    if($$params{projectChanged}){
        $$params{browseDir} = $project{directory};
        $$params{browseDirChanged} = 1;
    }    
    $$params{browseDir} or $$params{browseDir} = $startDir;
    if($$params{browseDirChanged} and $$params{browseDirItem}){
        if($$params{browseDirItem} eq '..'){
            $$params{browseDir} =~ s|(.+)/\w+$|$1|;
        } else {
            $$params{browseDir} .= '/'.$$params{browseDirItem};
        }
    }  
    push @masterInputs, (
        tablewrap(
            [$cgi->textarea(-name=>"browseDir",
                             -id=>"browseDir",
                             -title=>"Current working directory",
	                         -class=>"browseDir",
	                         -rows=>2,
                             -onchange=>'changeBrowseDir()'),
             "directory"],
        ),     
        browsePopups()     
    );
}
sub browsePopups {
    if( !defined($$params{browseDirList}) or $$params{browseDirChanged} ){
        my $ls = "ls -1 --file-type $$params{browseDir}";
        my $ssh = getSSH($ls);
        parseDirContents($ssh, \my@dirs, \my@files);
        @dirs = sort {lc($a) cmp lc($b)} @dirs;
        $$params{browseDirList} = join(",", @dirs);
        $$params{browseDirItem} = "";
        @files = sort {lc($a) cmp lc($b)} @files;
        $$params{browseFileList} = join(",", @files);
        $$params{browseDirChanged} and $$params{browseFileItem} = '';  
    }
    return (
        tablewrap(
            [$cgi->popup_menu(-name=>"browseDirItem",
                                -id=>"browseDirItem",
                                -values=>['', '..', split(",", $$params{browseDirList})],
                                -title=>"Move to a directory",
                                -onChange=>"updateField('browseDirChanged')",
                                -class=>"standard").
             $cgi->hidden(-name=>"browseDirList", -id=>"browseDirList").
             "\n<input type=hidden name=browseDirChanged id=browseDirChanged>\n",
             "cd"],
        ),      
        tablewrap(
            [$cgi->popup_menu(-name=>"browseFileItem",
                                -id=>"browseFileItem",
                                -values=>['', split(",", $$params{browseFileList})],
                                -title=>"Select a file to view or edit",
                                -onChange=>"updateField('browseFileChanged')",
                                -class=>"standard").
             $cgi->hidden(-name=>"browseFileList", -id=>"browseFileList").
             "\n<input type=hidden name=browseFileChanged id=browseFileChanged>\n",
             "file"],
        ),    
    );                  
}
sub parseDirContents {
    my ($ssh, $dirs, $files) = @_;
    foreach my $item(split("\n", $ssh)){
        $item =~ m|\@$| and next;  # ignore links
        if($item =~ m|(.+)\/$|){
            push @$dirs, $1;
        } else {
            push @$files, $item;
        }
    }
}
#-----------------------------------------------------------
sub browseTaskButtons {
    unless($$params{browseDir} and $$params{browseFileItem}){ 
        $$params{currentTask} = ""; 
        $$params{currentSubTask} = ""; 
        return;
    }
    $isSSH or return;
    my $isImage = isImageFile();   
    if($isImage){
        $$params{isImage}++;
    } else {
        setBrowseCurrentTask(\@fileTasks, \%fileTasks, 'head');
    }
}
sub setBrowseCurrentTask {
    my ($tasks, $taskHash, $defaultTask) = @_;
    $$params{currentTask} or $$params{currentTask} = $defaultTask;
    $$taskHash{$$params{currentTask}} or $$params{currentTask} = $defaultTask;
    addBrowseTaskButtons($tasks, $taskHash);   
}
sub addBrowseTaskButtons {
    my ($tasks, $taskHash) = @_;
    foreach my $task(@$tasks){
        push @taskButtons, taskButton($task, $task, $$taskHash{$task});
    }
}
#-----------------------------------------------------------
sub getBrowse {
    $$params{browseDir} or return;
    $$params{browseFileItem} or return;
    my $remoteFile = "$$params{browseDir}/$$params{browseFileItem}";
    my $imageExtension = isImageFile();   
    if($imageExtension){
        my $localImage = "user/$config{user}.$imageExtension";
        my $localFile = "$imagesDir/$localImage";
        if($isServer){
            qx|cp $remoteFile $localFile|;
        } else {
            scpDownload($localFile, $remoteFile);        
        }
        my $rand = int(rand(1E9));  # force re-acquisition of image file with same name
        @taskResults = "<img width=100% src=\"../q_images/$localImage?rand=$rand\">";
    } elsif($$params{currentTask} eq 'Edit'){
        saveFile('browse');
    } else {
        my $cmd = "$$params{currentTask} $remoteFile"; 
        my $displayCmd = $cmd;     
        $isServer or $cmd = '"'.$cmd.'"';
        my $results = "EXECUTING:\n$displayCmd\n\nRESULTS:\n".getSSH("$cmd 2>&1");
        open my $inH, "<", \$results;
        return parseQReturn($inH);           
    }
}
sub isImageFile {
    $$params{browseFileItem} or return;
    $$params{browseFileItem} =~ m|.+\.(.+)$| or return;
    return $imageExtensions{$1}
}
############################################################

1;

