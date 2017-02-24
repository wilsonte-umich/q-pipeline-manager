#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($params @taskResults @errorResults @taskInputs %project 
            $mobileHR $qURL $stylesheet $javascript %css %js $libDir);
my $mobileScript = "$libDir/mobile.pl";
my $isMobileScript = (-e $mobileScript);
$isMobileScript and require $mobileScript;
############################################################

############################################################
# common html functions
#-----------------------------------------------------------
sub tablewrap {      # wraps a series of elements into a single table row
    my (@refs) = @_; # each element provided as [input (required), label (optional), hidden (optional)]
    my $row = "<tr>\n";
    foreach my $ref(@refs){
        my ($input, $label, $hidden) = @$ref;
        $label or $label = "";
        $hidden or $hidden = "";
        $row .= "<td>$label&nbsp</td><td>$input$hidden</td>\n";
    }
    $row .= "</tr>\n";
}
sub tablewrap2 {     # wraps a series of elements into a single table row
    my (@refs) = @_; 
    my $row = "<tr>\n";
    foreach my $ref(@refs){
        foreach my $cell(@$ref){
            $cell or $cell = "";
            $row .= "<td>$cell</td>\n";
        }
    }
    $row .= "</tr>\n";
}
sub noOptions {
    push @taskInputs, "no options";
}
############################################################

############################################################
# job option parsing for submission to q
#-----------------------------------------------------------
sub checkRequiredOptions {
    my ($taskResults, @optionNames) = @_;
    foreach my $optionName(@optionNames){ 
        $$params{$optionName} or push @errorResults, "<p class=errorMessage>Error:  option '$optionName' is required<br><br>"; 
    }
}
sub setBooleanOptions {
    my ($options, @optionNames) = @_;
    foreach my $optionName(@optionNames){ $$params{$optionName} and push @$options, "--$optionName" }
}
sub setValueOptions {
    my ($options, @optionNames) = @_;
    foreach my $optionName(@optionNames){ $$params{$optionName} and push @$options, "--$optionName $$params{$optionName}" }
}
############################################################

############################################################
# parse q return lines
#-----------------------------------------------------------
sub parseQReturn {
    my ($qH, $parseSub) = @_;
    $qH or return "failed to retrieve q handle";
    push @taskResults, "
    <div id=statusDiv>
    <table>\n";    
    my $qType;
    while(my $line = <$qH>){ push @taskResults, taskRowWrap($line, $., $parseSub) }
    push @taskResults, "
    </table>
    </div>\n";
    close $qH; 
}
sub taskRowWrap {
    my ($cell, $rowID, $parseSub) = @_;
    $parseSub and &$parseSub($cell);
    maskTaskRow(\$cell);   
    my $nbsp = '&nbsp;';
    "<tr id=$rowID name=$rowID onmouseover=changeRowColor(this,true,'$rowID') onmouseout=changeRowColor(this,false,'$rowID')>".
    "<td class=status nowrap=nowrap>$cell<input type=hidden id=isClicked$rowID value=0 /></td></tr>\n";  
}
sub maskTaskRow{
    my ($cell) = @_;
    chomp $$cell;
    $$cell =~ s|~||g;
    $$cell =~ s|$project{directory}/|\.\./|g;
    $$cell =~ s| |&nbsp;|g;
    return 1;  
}
############################################################

############################################################
# general utilities
#-----------------------------------------------------------
sub resetPage {  # clear all form inputs
    foreach my $key(%$params){ $key and defined $$params{$key} and delete $$params{$key} }
    resetCache();
}
sub slurpFile {  # read the entire contents of a disk file into memory
    my ($file) = @_;
    local $/ = undef; 
    -e $file or return "";
    open my $inH, "<", $file or die "\ncould not open $file for reading: $!\n";
    my $contents = <$inH>; 
    close $inH;
    return $contents;
}
sub stripComments {  # just the important content of a file
    my ($lines) = @_;
    my @stripped;
    open my $inH, "<", $lines;
    while(my $line = <$inH>){
        chomp $line;
        $line =~ s/\r//g;  # dos2Unix
        $line =~ m/^\s*([^#]*)/; # strip leading white space and ignore comments
        $1 or next; # ignore lines with no content
        $line = $1;
        $line =~ s/\s+$//; # strip trailing white space
        push @stripped, "$line\n";
    }
    close $inH; 
    return join("", @stripped);
}
#-----------------------------------------------------------
sub test {  # used during development to understand how information is passed
    my ($request) = @_;
    my @html;
    push @html,
    qw(<!DOCTYPE html>
    <html>
    <head>
    <title>q Remote</title>
    <table>
    test<br>);
    push @html, "<tr><td>-----------<br></td></tr><tr><td>params<br></td></tr>";
    foreach my $key (keys %$params){ push @html, "<tr><td>$key</td><td>$$params{$key}</td></tr>" }
    push @html, "<tr><td>-----------<br></td></tr><tr><td>request<br></td></tr>";
    foreach my $key (keys %$request){ push @html, "<tr><td>$key</td><td>$$request{$key}</td></tr>" }
    push @html, "<tr><td>-----------<br></td></tr><tr><td>request_headers<br></td></tr>";
    foreach my $key (keys %{$$request{_headers}}){ push @html, "<tr><td>$key</td><td>$$request{_headers}{$key}</td></tr>" }
    push @html,$request->url,
    qw(</table>
    </body>
    </html>);
    return join("\n", @html);
}
############################################################


############################################################
# set screen formatting based on client device type
#-----------------------------------------------------------
sub isMobile {
    return $isMobileScript ? checkMobile() : 0;
}
#-----------------------------------------------------------
sub setScreenFormatting {
    my ($isMobile) = @_;
    if($isMobile){
        $stylesheet = $css{common} . $css{mobile}; 
        $javascript = $js{common} . $js{mobile};
        $mobileHR = "<hr>"; 
        $qURL = "http://tewlab.org/q/mobile.html";
    } else {
        $stylesheet = $css{common} . $css{desktop}; 
        $javascript = $js{common} . $js{desktop}; 
        $mobileHR = ""; 
        $qURL = "http://tewlab.org/q";
    }
}
############################################################

1;

