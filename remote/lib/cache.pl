#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($cgi $params @taskResults $isServer $debug);
our (@cacheButtons, $isCacheRequest); 
my ($currentIndex, $maxCache, @cache, %lastPageState) = (0, 10);
my %pageCacheParams = (
    queue => [qw(project masterClass masterClassList masterName masterNameList isLocked currentTask currentSubTask)],
    shell => [qw(project browseDir browseDirList browseFileList browseFileItem currentTask shellCommand)]
);
############################################################

############################################################
# manage page cache
#-----------------------------------------------------------
sub moveToCache {  # if requested move to cached page
    $isCacheRequest = undef;
    $$params{cacheDirection} or return;  # not a cache request
    $isCacheRequest = 1;
    $currentIndex += $$params{cacheDirection};
    my ($request, $taskResults) = @{$cache[$currentIndex]};
    @taskResults = @$taskResults;
    $cgi = CGI->new($request->content);
    $params = $cgi->Vars;
    $$params{executing} = undef;  # never re-execute on a cache call; could result in re-queuing jobs
}
sub setCache {
    my ($request) = @_;
    $isServer and return;  # caching only meaningful for ssh-driven daemon
    unless($isCacheRequest){  # just read from cache, no new page was generated
        $currentIndex == 0 or @cache = $cache[$currentIndex..$#cache];  # clear forwards when new page
        unshift @cache, [$request, [@taskResults]];  # store current page in cache
        $#cache > $maxCache and pop @cache;  # enforce cache limit    
        $currentIndex = 0;  
    }
    setCacheButtons();
}
sub setCacheButtons {
    my $value = $#cache;
    my ($activeButton, $silencedButton) = ('activeButton', 'silencedButton');
    my ($backButtonClass, $backButtonScript) = ($activeButton, "setCacheDirection('1')");
    my ($forwardButtonClass, $forwardButtonScript) = ($activeButton, "setCacheDirection('-1')");
    $currentIndex < $#cache or ($backButtonClass, $backButtonScript) = ($silencedButton, "");
    $currentIndex > 0 or ($forwardButtonClass, $forwardButtonScript) = ($silencedButton, "");
    @cacheButtons = (
    "<input type=hidden name=cacheDirection id=cacheDirection>",
    "<input type=button class=$backButtonClass value=Back ".
            "onClick=$backButtonScript title=\"Go to previous page\">&nbsp&nbsp", 
    "<input type=button class=$forwardButtonClass value=Forward ".
            "onClick=$forwardButtonScript title=\"Go to next page\">&nbsp&nbsp" 
    );
}
sub resetCache {
    ($currentIndex, $maxCache) = (0, 10);
    @cache = ();
}
#-----------------------------------------------------------
sub savePageState { # separate cache for handling Queue to Shell page type switch
    $$params{pageType} or return;  # this cache goes with page to client
    foreach my $param(@{$pageCacheParams{$$params{pageType}}}){ 
        my $value = $$params{$param};
        defined $value or next;
        push @cache, "$param=$value" 
    }
    $$params{"pageCache_$$params{pageType}"} = join(";", @cache);    
    $cgi->hidden(-name=>'pageCache_queue', -id=>'pageCache_queue').
    $cgi->hidden(-name=>'pageCache_shell', -id=>'pageCache_shell')
}
sub getPageState {
    $$params{pageTypeChanged} or return;
    $$params{pageType} or return;
    my $pageCache = $$params{"pageCache_$$params{pageType}"};
    $pageCache or return;
    my @cache = split(";", $pageCache);
    foreach my $cache(@cache){
        my ($param, $value) = split("=", $cache);
        $$params{$param} = $value;
    }
}
############################################################

1;

