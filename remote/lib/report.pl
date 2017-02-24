#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @taskInputs @taskResults);
############################################################

############################################################
# q status web interface
#-----------------------------------------------------------
sub getReport {
    my @options = ("--job $$params{getReport}");
    my $qH = getQH(['report', @options], ['script', @options], ['environment', @options]); 
    $qH or return "failed to retrieve q handle";
    push @taskResults, 
    "<div id=statusDiv>
    <table>\n";
    while(my $line = <$qH>){       
        push @taskResults, reportRowWrap(\$line);
    }
    push @taskResults, "
    </table>
    </div>\n";
    close $qH;  
}
sub reportRowWrap {
    my ($line) = @_;
    maskTaskRow($line) or return "";      
    highlightReportLine($line);    
    "<tr><td class=status nowrap=nowrap>$$line</td></tr>\n";     
}
sub highlightReportLine {
    my ($line) = @_;
    $$line =~ m|\:&nbsp;job\:&nbsp;\d+| and $$line = "<b>$$line</b>";   
    $$line =~ s|(\$\w+)|<font color=maroon>$1</font>|g;
    $$line =~ s|(#.*)|<font color=darkgreen>$1</font>|g;  
}
############################################################

1;

