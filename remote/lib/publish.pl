#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @taskInputs @taskResults $imagesDir $isServer %config);
our $zipFile = undef;
############################################################

############################################################
# publish pipeline to static HTML
#-----------------------------------------------------------
sub publishInputs {
    push @taskInputs, (
        $cgi->checkbox(-name=>'long',
                       -title=>"Include jobID and start time information when publishing"),
        "<div class=vSpacer></div>
        <table>",
         tablewrap(
            [$cgi->textfield(-name=>"title",
                             -id=>"title",
                             -title=>"Title of the published html report (overridden by 'publishTitle' instructions)",
	                         -class=>"wide"),
             "title"],
        ),   
        tablewrap(
            [$cgi->textfield(-name=>"intro-file",
                             -id=>"intro-file",
                             -title=>"HTML overview file to include when publishing (overridden by 'publishIntro' instructions)",
	                         -class=>"wide"),
             "intro-file"],
        ),       
        tablewrap(
            [$cgi->textfield(-name=>"mask",
                             -id=>"mask",
                             -title=>"Comma-delimited list of strings to mask when publishing",
	                         -class=>"wide"),
             "mask"],   
        ),  
        tablewrap(
            [$cgi->textfield(-name=>"out-dir",
                             -id=>"out-dir",
                             -title=>"Publish output directory (overridden by 'publishDir' instructions)",
	                         -class=>"wide"),
             "out-dir"], 
        ),           
        "</table>",
    );
} 
sub getPublish {
    setBooleanOptions(\my@options, qw(long));
    setValueOptions(\@options, qw(mask title intro-file out-dir));
    my $qH = getQH(['publish', @options]); 
    parseQReturn($qH, \&parsePublishZip); 
    if($zipFile){
        $zipFile =~ m|.+/(.+\.tar\.gz)| or return;
        my $fileName = $1;
        my $userPath = "user/$config{user}";
        my $localDir = "$imagesDir/$userPath";
        mkdir $localDir;
        my $localPath = "$userPath/$fileName";
        my $localFile = "$imagesDir/$localPath";
        if($isServer){
            qx|cp $zipFile $localFile|;  
        } else {
            scpDownload($localFile, $zipFile);            
        }
        my $rand = int(rand(1E9));
        push @taskResults, 
            "<div id=statusDiv>
             <table>
             <tr><td class=status nowrap=nowrap><a href=\"../q_images/$localPath?rand=$rand\")>get zipped HTML report</a></td></tr>
             </table>
             </div>\n";              
    }
}
sub parsePublishZip {
    my ($line) = @_;
    $line =~ m|(\S+\.tar\.gz)| and $zipFile = $1;
}
############################################################

1;

