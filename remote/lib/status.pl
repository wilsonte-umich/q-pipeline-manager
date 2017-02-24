#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @taskInputs @taskResults $debug);
my @statusColumns = qw(user job_name job_ID exit_status start_time wall_time maxvmem);
my @sortOrders = qw(asc desc);
my @filterOperators = qw(= != ~ !~ > <);
############################################################

############################################################
# q status web interface
#-----------------------------------------------------------
sub statusInputs {
    push @taskInputs, (
        $cgi->checkbox(-name=>'archive',
                       -title=>"Show the most recent archive instead of the current status"),
        "<div class=vSpacer></div>
        <table>",
        tablewrap(
            [$cgi->popup_menu(-name=>"sortColumn", 
                             -values=>["", @statusColumns],
                             -title=>"Column by which status lines should be sorted",
                             -class=>"column"),
             "sort"],
            [$cgi->popup_menu(-name=>"sortOrder", 
                             -values=>[@sortOrders],
                             -title=>"Sort order (ascending or descending)",
                             -class=>"operator")]
        ),
        tablewrap(
            [$cgi->popup_menu(-name=>"filterColumn", 
                             -values=>["", @statusColumns],
                             -title=>"Column by which status lines should be filtered",
                             -class=>"column"),
             "filter"],
            [$cgi->popup_menu(-name=>"filterOperator", 
                             -values=>[@filterOperators],
                             -title=>"Filter operator",
                             -class=>"operator")],
            [$cgi->textfield(-name=>"filterValue",
                            -title=>"Filter status lines by this value",
	                        -class=>"standard")]
        ),           
        $cgi->hidden(-name=>"job", -id=>"job"),
        "</table>",
    );
}
sub getStatus {
    setBooleanOptions(\my@options, qw(archive));
    $$params{sortColumn} and push @options, "--sort $$params{sortColumn}:$$params{sortOrder}";
    $$params{filterColumn} and $$params{filterValue} and 
        push @options, "--filter '$$params{filterColumn}$$params{filterOperator}$$params{filterValue}'";
    my $qH = getQH(['status', @options]); 
    $qH or return "failed to retrieve q handle";
    push @taskResults, "
    <div id=statusDiv>
    <table>\n";
    my $qType;
    while(my $line = <$qH>){
        my @line = split(/\s+/, $line);        
        $line =~ m|^submitted\s+| and $qType = $line[1];
        my $jobID;
        $line[2] and $line[2] =~ m|^(\d+)$| and $jobID = $1;
        !$jobID and $line[3] and $line[3] =~ m|^(\d+)$| and $jobID = $1;
        $jobID or $jobID = "deadRow";
        push @taskResults, statusRowWrap($line, $jobID, 'processStatusListClick');
    }
    push @taskResults, "
    </table>
    </div>\n";
    close $qH;  
    $$params{isStatus} = 1;
}
sub statusRowWrap {
    my ($cell, $rowID) = @_;
    maskTaskRow(\$cell) or return "";   
    my $nbsp = '&nbsp;';
    my $link = "<a href=javascript:getReport('$rowID')>$rowID</a>";
    $cell =~ s|$nbsp$rowID$nbsp|$nbsp$link$nbsp|;
    "<tr id=$rowID name=$rowID onmouseover=changeRowColor(this,true,'$rowID') onmouseout=changeRowColor(this,false,'$rowID')>".
    "<td class=status nowrap=nowrap onclick=processStatusListClick(event,this,'$rowID')>".
    "$cell<input type=hidden id=isClicked$rowID value=0 /></td></tr>\n"; 
}
############################################################

1;

