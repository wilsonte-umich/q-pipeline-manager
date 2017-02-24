#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params %project @taskResults @errorResults
            %isLocked %template %help $isServer $debug);
our @masterInputs;
our $updateMasterClasses = undef;
our $updateMasterNames = undef;
############################################################

############################################################
# master selection inputs
#-----------------------------------------------------------
sub masterParameters {  
    @masterInputs = ();
    $$params{pageType} eq 'queue' or return;
    $project{directory} or return;
    push @masterInputs,
    tablewrap2(['master class&nbsp;',
                 masterClassPopup().
                 "\n<input type=hidden name=masterClassChanged id=masterClassChanged>\n",
                 "&nbsp;&nbsp;<a href=javascript:newMaster('class') title=\"Create a new master class directory on server\">New</a>".
                 "\n<input type=hidden name=newMasterclass id=newMasterclass>\n"]);
    $$params{masterClass} or return;
    push @masterInputs,
    tablewrap2(['master name&nbsp;',
                 masterNamePopup()."\n<input type=hidden name=masterNameChanged id=masterNameChanged>\n",
                 "&nbsp;&nbsp;<a href=javascript:newMaster('') title=\"Create a new master file on server\">New</a>".
                 "\n<input type=hidden name=newMaster id=newMaster>\n"]);       
}
#-----------------------------------------------------------
sub masterClassPopup {
    $$params{masterClassChanged} and $$params{masterName} = "";
    if(!$$params{masterClassList} or ($$params{projectChanged} or $updateMasterClasses)){
        my $ssh;
        if($isServer){
            foreach my $dir(<$project{directory}/masters/*>){
                -d $dir or next;
                $ssh .= "$dir/\n";
            }
        } else {
            my $ls = "ls -1 -d $project{directory}/masters/*/";
            $ssh = getSSH($ls);
        }
        my @masterClasses = map { maskMasterClass($_) } split("\n", $ssh);  
        @masterClasses = sort {lc($a) cmp lc($b)} @masterClasses;
        $$params{masterClassList} = join(",", @masterClasses);
        $updateMasterClasses = undef;
        $updateMasterNames = 1;
        $$params{masterName} = "";
    }
    return $cgi->popup_menu(-name=>"masterClass",
                            -id=>"masterClass",
                            -values=>["", split(",", $$params{masterClassList})],
                            -title=>"The class of master to create, run or monitor",
                            -onChange=>"updateField('masterClassChanged')",
                            -class=>"standard").
           $cgi->hidden(-name=>"masterClassList", 
                        -id=>"masterClassList");
}
sub maskMasterClass{
    my ($masterClass) = @_;
    $masterClass =~ s|^$project{directory}/masters/||;
    $masterClass =~ s|/.*$||;
    return $masterClass;
}
#-----------------------------------------------------------
sub masterNamePopup {
    if(!$$params{masterNameList} or ($$params{masterClassChanged} or $updateMasterNames)){   
        my $ssh;
        if($isServer){
            foreach my $dir(glob("$project{directory}/masters/$$params{masterClass}/*.q*")){
                -d $dir and next;
                $ssh .= "$dir/\n";
            }
        } else {
            my $ls = "ls -1 -f $project{directory}/masters/$$params{masterClass}/*.q*";
            $ssh = getSSH($ls);
        }
        my @masterNames = ();
        $$params{isLocked} = "";
        foreach my $masterName(split("\n", $ssh)){
            $masterName = maskMasterName($masterName);
            $masterName and push @masterNames, $masterName;
        } 
        @masterNames = sort {lc($a) cmp lc($b)} @masterNames;
        $$params{masterNameList} = join(",", @masterNames);
        $updateMasterNames = undef;
    }
    return $cgi->popup_menu(-name=>"masterName",
                             -id=>"masterName", 
                             -values=>["", split(",", $$params{masterNameList})],
                             -title=>"The name of the master to create, run or monitor",
                             -onChange=>"updateField('masterNameChanged')",
                             -class=>"standard").
           $cgi->hidden(-name=>"masterNameList", 
                        -id=>"masterNameList").
           $cgi->hidden(-name=>"isLocked", 
                        -id=>"isLocked");
}
sub maskMasterName{
    my ($masterName) = @_;
    my $masterFile = $masterName;
    $masterName =~ s|^$project{directory}/masters/$$params{masterClass}/||;
    $masterName =~ s|/s*$||;
    if($masterName =~ m|(.+)\.lock|){
        $$params{isLocked} .= "$1,";
    } elsif($masterName =~ m|(.+)\.q\.template|){
        $template{$masterName} = $masterFile;
    } elsif($masterName =~ m|(.+)\.q\.help|){
        $help{$masterName} = $masterFile;
    } else {
        return $masterName;
    }
    return undef;
}
############################################################

############################################################
# create new master classes and masters
#-----------------------------------------------------------
sub newMasterClass {
    $$params{pageType} eq 'queue' or return;
    $$params{newMasterclass} or return;
    $project{directory} or return;
    $$params{masterClassList} or return;
    my %masterClasses = map { $_ => 1 } split(",", $$params{masterClassList});
    if( $masterClasses{$$params{newMasterclass}} ){
        @errorResults = ("<p class=errorMessage>Error: master class $$params{newMasterclass} already exists");
    } else {
        my $newMasterDir = "$project{directory}/masters/$$params{newMasterclass}";
        my $mkdir = "mkdir $newMasterDir ; chmod -R ug+rw $newMasterDir";
        $isServer or $mkdir = "'$mkdir'";
        @errorResults = (getSSH($mkdir));
    }
    @errorResults and return;
    $$params{projectChanged} = 1;
    $$params{masterClassChanged} = 1;
    $updateMasterClasses = 1;
    $$params{masterClass} = $$params{newMasterclass};
}
sub newMasterName {
    $$params{pageType} eq 'queue' or return;
    $$params{newMaster} or return;
    $project{directory} or return;
    $$params{masterClass} or return;
    $$params{masterNameList} or $$params{masterNameList} = "";   
    my %masterNames = map { $_ => 1 } split(",", $$params{masterNameList});
    $$params{newMaster} =~ m|\.q$| or $$params{newMaster} = "$$params{newMaster}.q"; 
    if( $masterNames{$$params{newMaster}} ){
        @errorResults = ("<p class=errorMessage>Error: master $$params{newMaster} already exists");
    } else {
        my $newMasterFile = "$project{directory}/masters/$$params{masterClass}/$$params{newMaster}";    
        my $templateFile = getFile('template');
        my $createFile = "touch $newMasterFile ; cp $templateFile $newMasterFile 2>/dev/null ; chmod ug+rw $newMasterFile";
        $isServer or $createFile = "'$createFile'";
        @errorResults = (getSSH($createFile));
    }
    @errorResults and return; 
    $updateMasterNames = 1;
    $$params{masterName} = $$params{newMaster};
    delete $$params{currentTask};
    delete $$params{currentSubTask};
}
############################################################

1;

