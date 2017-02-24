#!/usr/bin/perl
package q_remote;
use strict;
use warnings;
use Cwd (qw(abs_path));
use HTTP::Daemon;
use CGI;  

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw(%config);
our $scriptName = "daemon.pl";
our $script = abs_path($0);
($script =~ m|(.+)/$scriptName| or $script =~ m|(.+)\\$scriptName|) or die "error recovering script directory\n";
our $scriptDir = $1;
our $libDir = "$scriptDir/lib";
our $imagesDir = "$scriptDir/images";
our $transferDir = "$scriptDir/transfer";
############################################################

############################################################
# read command line options
#-----------------------------------------------------------
my %optionInfo = (# [shortOption, isValueOption]          
    'no-agent'=> ["n", undef],
    'config'=>   ["c", 1]
);
my %longOptions = map { ${$optionInfo{$_}}[0] => $_ } keys %optionInfo; 
our %options;
while (my $optionList = shift @ARGV){  
    if($optionList =~ m/^\-\-(.+)/){ # long option formatted request
        my $longOption = $1;
        defined $optionInfo{$longOption} or die "'$longOption' is not a recognized option\n"; 
        setOption($longOption);
    } elsif ($optionList =~ m/^\-(.+)/){ # short option formatted request
        foreach my $shortOption(split('', $1)){
            my $longOption = $longOptions{$shortOption};
            defined $longOption or die "'$shortOption' is not a recognized option\n"; 
            setOption($longOption);
        }
    } else {
        die "malformed option list\n"; 
    }
}
$options{config} or $options{config} = 'daemon.conf';
#-----------------------------------------------------------
sub setOption { # check and set option request                
    my ($longOption) = @_;
    my $value = ${$optionInfo{$longOption}}[1] ? shift @ARGV : 1;
    defined $value or die "missing value for option '$longOption'\n";
    $value =~ m/^\-/ and die "missing value for option '$longOption'\n";
    $options{$longOption} = $value;  
}
############################################################

############################################################
# load required files
#-----------------------------------------------------------
require "$libDir/config.pl";
require "$libDir/ssh.pl";
require "$libDir/cache.pl";
require "$libDir/utilities.pl";
require "$libDir/main.pl";
require "$libDir/project.pl";
require "$libDir/master.pl";
require "$libDir/task.pl";
#-----------------------------------------------------------
require "$libDir/submit.pl";
require "$libDir/status.pl";
require "$libDir/clear.pl";
require "$libDir/delete.pl";
require "$libDir/resubmit.pl";
require "$libDir/report.pl";
require "$libDir/lock.pl";
require "$libDir/purge.pl";
require "$libDir/protect.pl";
require "$libDir/backup.pl";
require "$libDir/edit.pl";
require "$libDir/publish.pl";
#-----------------------------------------------------------
require "$libDir/queue.pl";
require "$libDir/browse.pl";
require "$libDir/shell.pl";
#-----------------------------------------------------------
loadConfig(); 
loadProjects();
############################################################

############################################################
# start the various parent and child processes
#-----------------------------------------------------------
my $url = "http://localhost:$config{localPort}/q/daemon.pl";
my $isAgent = ($config{agentCommand} and !$options{'no-agent'});
my $isPageant = ($config{agentCommand} and lc($config{agentCommand}) =~ m|pageant|);
$isAgent and ($isPageant or startSshAgent($isPageant));  # ssh-add returns, so is inline
my $pid = fork();
if($pid or !(defined $pid)){
    startDaemon();
} else {
    $isAgent and ($isPageant and startSshAgent($isPageant)); # Pageant does not return
    startBrowser();  # won't start if $isAgent and $isPageant
}
############################################################

############################################################
# start the q remote http daemon on localhost
#-----------------------------------------------------------
sub startDaemon {
    my $daemon = HTTP::Daemon->new(LocalAddr=>'127.0.0.1', # ReuseAddr => 1, ReusePort
                                   LocalPort=>$config{localPort})
                 or die "could not start web daemon on $config{localPort}\n";
    print "\nPlease point your web browser to $url\n\nCtrl-C to exit Daemon\n";                          
    while (my $client = $daemon->accept) {
        while (my $request = $client->get_request) {
            my $host = lc($$request{_headers}{host});   
            my ($page) = $request->url->path =~ m|.*/q/(.+)|;
            my $isImage; 
            unless ($page){
                ($page) = $request->url->path =~ m|.*/q_images/(.+)|; 
                $page and $isImage = 1;
            }
            if(!($host =~ m|^localhost| or $host =~ m|^127.0.0.1|)){       
                $client->send_error();  
            } elsif (!$page) {  
                $client->send_error();
            } elsif ($page eq 'daemon.pl') {
                my $html = q_remote($request);
                my $response = HTTP::Response->new('200', 'OK', [ 'Content-Type' => 'text/html'], $html);   
                $client->send_response( $response );
            } elsif ($isImage) {
                $client->send_file_response( "$imagesDir/$page" );
            } else {
                $client->send_error();
            }
        }
        $client->close();
        undef($client);
    }
}
############################################################

############################################################
# use ssh key agent to prompt for passphrases
#-----------------------------------------------------------
sub startSshAgent {
    my ($isPageant) = @_;
    print "\nloading keyfiles into ssh agent   $config{agentCommand}\n";
    if($isPageant){
        $config{keyFile} or die "config error: keyFile required when using Pageant\n";
        exec($config{agentCommand}, $config{keyFile});
        exit;
    } else {
        $config{keyFile} or $config{keyFile} = "";
        system("$config{agentCommand} $config{keyFile}");   
    }
}
############################################################

############################################################
# open the required browser window
#-----------------------------------------------------------
sub startBrowser {
    $config{browserCommand} or return;
    sleep 1; # give the daemon a moment to fire up
    print "\nopening browser window\n";
    exec($config{browserCommand}, $url);
}
############################################################


