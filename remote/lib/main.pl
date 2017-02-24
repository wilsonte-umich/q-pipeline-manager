#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($scriptDir $libDir @cacheButtons %project @projectInputs @masterInputs 
            @taskButtons @subTaskButtons @taskInputs @taskResults 
            $sshError $isServer $debug);
our ($cgi, $params, @refreshButtons, @pageButtons, @goButtons, 
    $isMobile, $stylesheet, $javascript, $mobileHR, $qURL);
our @errorResults = ();
our $resetting = undef;
our %css = map { $_ => slurpFile("$libDir/q_remote_$_.css") } qw(common desktop mobile);
our %js = map { $_ => slurpFile("$libDir/q_remote_$_.js") } qw(common desktop mobile);
my %isDryRun = map { $_ => 1 } qw(submit extend resubmit clear delete rollback purge 
                                  protect unprotect backup restore master help template);
my %isEdit = map { $_ => 1 } qw(master help template);
my $version = slurpFile("$scriptDir/VERSION");
chomp $version;
############################################################

############################################################
# assemble the parts of the q remote web page
#-----------------------------------------------------------
sub q_remote {
    my ($request) = @_;   
    $isMobile = isMobile();
    setScreenFormatting($isMobile);
    $cgi = $request ? CGI->new($request->content) : CGI->new;
    $params = $cgi->Vars; 
    @errorResults = ();
    $isServer or moveToCache();    
    $$params{reset} and resetPage();  
    $$params{pageType} or $$params{pageType} = 'queue';
    getPageState();
    setProject();    
    checkSSH(); 
    setIsLocked(); 
    taskResults();
    newMasterClass();
    newMasterName();
    toggleLock();    
    refreshButtons();  # order is important, can't commit to client until all are done
    masterParameters();
    setIsLocked();  # yes, do it twice   
    taskButtons();
    subTaskButtons(); 
    taskInputs();
    pageButtons();
    goButtons();
    setCache($request);     
    @errorResults and @taskResults = ();     
    return generateHtml();
}
############################################################

############################################################
# generate the final html
#-----------------------------------------------------------
sub generateHtml {
    my @html;
    push @html,
    "<!DOCTYPE html>
    <html>
    <head>
    <title>Remote</title>
    <meta name=description content=q_remote />
    <meta name=author content=\"Thomas E. Wilson\" />
    <meta name=copyright content=\"2011© Thomas E. Wilson, University of Michigan\" />
    <meta http-equiv=content-type content=\"text/html;charset=UTF-8\" />
    <link rel=icon type=image/png href=../q_images/q_icon_50.png>
    <style type=text/css>
    $stylesheet
    </style>
    <script type=text/javascript>
    $javascript
    </script>
    </head>
    <body onload=setPage() onresize=setSizes(0)>
    <form name=q_remote method=post>
    <div id=inputs class=inputs>
    <div class=icon><a title=\"$version (click for extended help)\" target=\"q_frame\" href=\"$qURL\"><img src=../q_images/q_icon_35.png></a></div> 
    @refreshButtons
    @cacheButtons 
    <div class=vSpacer></div>
    @pageButtons
    <hr>
    <table> 
    @projectInputs
    @masterInputs
    </table>\n",
    $sshError ? "<p class=errorMessage>$sshError</p>" : "",    
    "@taskButtons
    @subTaskButtons
    @taskInputs
    @goButtons
    </div>
    $mobileHR
    <div id=sshReturn class=sshReturn>
    @errorResults
    @taskResults
    </div>
    </form>
    </body>
    </html>";
    return join("\n", @html);
}
############################################################

############################################################
# generate the page-level action buttons
#-----------------------------------------------------------
sub refreshButtons {
    @refreshButtons = (
        "<input type=submit class=activeButton value=Refresh title=\"Reload the current page\">&nbsp&nbsp",
        "<input type=button class=activeButton value=Reset onClick=resetPage() ".
                "title=\"Forget all settings and start fresh\">&nbsp&nbsp",
        "<input type=hidden name=reset id=reset>"
    );
}
#-----------------------------------------------------------
sub pageButtons {
    @pageButtons = (
        pageButton('queue', 'Queue', 'Use q to access, execute and manage pipelines'),
        pageButton('shell', 'Shell', 'Browse directories, edit files and execute commands on server'),
        $cgi->hidden(-name=>'pageType', -id=>'pageType'),
        "<input type=hidden name=pageTypeChanged id=pageTypeChanged>",
        savePageState()
    );
}
sub pageButton {
    my ($pageType, $label, $title) = @_;
    return "<input type=button class=activeButton id=$pageType"."Button value=\"$label\"
                   onClick=setPageType('$pageType') title=\"$title\">&nbsp&nbsp",
}
#-----------------------------------------------------------
sub goButtons {
    @goButtons = ();
    @taskInputs or return;
    push @goButtons, (
        "<hr>
        <input type=hidden name=executing id=executing >
        <input type=hidden name=isDryRun id=isDryRun >"
    );
    my ($goLabel, $goHelp, $dryRunLabel, $dryRunHelp);
    my $isEdit =  ($$params{currentTask} eq '_edit' or $$params{currentTask} eq 'Edit') ? 
                  ($isEdit{$$params{currentSubTask}} or $$params{currentTask} eq 'Edit') : 
                  undef;
    if($$params{isImage}){               
        $goLabel = 'Display';  
        $goHelp = 'Show the file image';           
    } elsif($isEdit){
        $goLabel = 'Save';  
        $goHelp = 'Save the edited file to host';  
        $dryRunLabel = 'Cancel';   
        $dryRunHelp = 'Discard any file edits'; 
    } elsif($$params{pageType} eq 'shell'){
        $goLabel = 'Go!';  
        $goHelp = 'Send command to host for execution';  
    } else {
        $goLabel = 'Go!';  
        $goHelp = 'Send task to host for execution by q';  
        $dryRunLabel = 'Dry Run';   
        $dryRunHelp = 'Check syntax and content of configured task'; 
    }  
    if($isEdit or $isDryRun{$$params{currentSubTask}}){
        push @goButtons, (
            "<input type=button class=activeButton id=dryRun value=\"$dryRunLabel\"
                    onClick=executeTask(1,'$$params{currentTask}','$$params{currentSubTask}') title=\"$dryRunHelp\">"
        )
    }
    push @goButtons, (
        "<input type=button class=activeButton id=execute value=\"$goLabel\"
                onClick=executeTask(0,'$$params{currentTask}','$$params{currentSubTask}') title=\"$goHelp\">"
    )
}
############################################################

1;

