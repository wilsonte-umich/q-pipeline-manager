#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @taskInputs @taskResults %project $scriptDir $isServer $debug);
my $singleQuote = "__SINGLEQUOTE__";
my $newLine = "__NEWLINE__";
############################################################

############################################################
# handle output file backups
#-----------------------------------------------------------
sub editInputs {
    my ($fileType) = @_;
    push @taskInputs, "<p class=input>Be sure to Save any file edits, or they will be lost.</p>";
    getFileContents($fileType);
}
sub getFileContents {
    my ($fileType) = @_;
    @taskResults = ();
    my $file = getFile($fileType);
    $file or return;
    my $fileContents = getSSH("cat $file");
    push @taskResults, (
        "<textarea name=editBox id=editBox class=editBox>$fileContents</textarea>",
        "<input type=hidden name=maskedEditBox id=maskedEditBox>"
    );
    $$params{editingFile}++;
}
sub saveFile {
    my ($fileType) = @_;
    $$params{isDryRun} and return;  # Cancel was clicked
    my $remoteFile = getFile($fileType);
    $remoteFile or return;
    my $localFile = "$scriptDir/transfer.tmp";    
    $$params{maskedEditBox} or $$params{maskedEditBox} = "";  
    $$params{maskedEditBox} =~ s/$singleQuote/'/g;
    $$params{maskedEditBox} =~ s/$newLine/\n/g;
    my $outFile = $isServer ? $remoteFile : $localFile;
    open my $outH, ">", $outFile or return;
    print $outH $$params{maskedEditBox};
    close $outH;
    $isServer or scpUpload($localFile, $remoteFile);
    getSSH("chmod ug+rw $remoteFile");
    $isServer or unlink $localFile;
}
sub getFile {
    my ($fileType) = @_;
    my %files;
    if($fileType eq 'browse'){
        $files{browse} = "$$params{browseDir}/$$params{browseFileItem}";
    } else {
        my $masterClass = "$project{directory}/masters/$$params{masterClass}";
        $files{master} = "$masterClass/$$params{masterName}";
        $files{template} = "$masterClass/$$params{masterClass}.q.template";
        $files{help} = "$masterClass/README";
    }
    return $files{$fileType};
}
############################################################

1;

