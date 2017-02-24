#!/usr/bin/perl
package q_remote;
use strict;
use warnings;
use Cwd (qw(abs_path));
use CGI;  

############################################################
# declare variables
#-----------------------------------------------------------
our $scriptName = "server.pl";
our $script = abs_path($0);
($script =~ m|(.+)/$scriptName| or $script =~ m|(.+)\\$scriptName|) or die "error recovering script directory\n";
our $scriptDir = $1;
our $libDir = "$scriptDir/lib";
our $imagesDir = "$scriptDir/images";
our $isServer = 1;
our $debug = "<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>";
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
# serve the page
#---------------------------------------------------------
print CGI::header();
print q_remote();
print "$debug\n";
############################################################

