#!/usr/bin/perl
use strict;
use warnings;
use Cwd(qw(abs_path));

#========================================================================
# 'configure.pl' sets up q for use on a host server system
# and generates the q program target
#========================================================================

#========================================================================
$| = 1;
#------------------------------------------------------------------------
# get the current version
#------------------------------------------------------------------------
my $script = abs_path($0);
$script =~ m|(.*)/configure.pl$| or die "fatal error: could not establish q program directory\n";
my $qDir = $1;
my $versionFile = "$qDir/VERSION";
my $versionError = "fatal error: could not establish q software VERSION\n";
open my $versionH, "<", $versionFile or die $versionError;
my $version = <$versionH>;
close $versionH;
($version) = $version =~ m|q\s+(\d+\.\d+\.\d+)|;
$version or die $versionError;
#------------------------------------------------------------------------
# check the license
#------------------------------------------------------------------------
my $licenseFile = "$qDir/LICENSE";
-e $licenseFile or die "fatal error: could not find q software LICENSE\n";
#------------------------------------------------------------------------
# get path to perl
#------------------------------------------------------------------------
print "reading system information ";
print ".";
my $perlPath = qx|which perl 2>/dev/null|;
$perlPath or die "fatal error\n".
                 "q requires that perl be available on the system\n";
chomp $perlPath;
#------------------------------------------------------------------------
# get path to bash
#------------------------------------------------------------------------
print ".";
my $bashPath = qx|which bash 2>/dev/null|;
$bashPath or die "fatal error\n".
                 "q requires that the bash shell be available on the system\n".
                 "http://www.gnu.org/software/bash/\n";
chomp $bashPath;        
#------------------------------------------------------------------------
# check path to /usr/bin/time
#------------------------------------------------------------------------
print ".";
my $timePath = qx|which /usr/bin/time 2>/dev/null|;
$timePath or die "fatal error\n".
                 "q requires that the GNU time utility (/usr/bin/time) be available on the system\n".
                 "http://www.gnu.org/software/time/\n";
chomp $timePath;  
#------------------------------------------------------------------------
# get time version and adjust memory correction accordingly
#------------------------------------------------------------------------
print ".";
my $timeError = "$timePath is not a valid installation of the GNU time utility\n";
my $timeVersion = qx/$timePath --version 2>&1 | head -n1/; 
chomp $timeVersion;
$timeVersion =~ m/GNU/ or die $timeError;
$timeVersion =~ m/time/ or die $timeError; 
my @tvf = split(/\s+/, $timeVersion);
$timeVersion = $tvf[$#tvf];
$timeVersion or die $timeError;  
$timeVersion eq "UNKNOWN" and $timeVersion = 1.8;
#($timeVersion and $timeVersion =~ s/GNU time //) or die "$timePath is not a valid installation of the GNU time utility\n";    
my $memoryCorrection = $timeVersion > 1.7 ? 1 : 4; # for time <= v1.7, account for the known bug that memory values are 4-times too large
my $memoryMessage = $timeVersion > 1.7 ? "" : "q: !! maxvmem value above is 4-fold too high due to known bug in GNU time utility !!"; 
#------------------------------------------------------------------------
# discover the job scheduler in use on the system
#------------------------------------------------------------------------
print ".";
my $qType = 0;  # no scheduler, will require submit option -e
if(isShellCommand('qhost')) {
    $qType = 'SGE';          
} elsif(isShellCommand('showq')) {
    $qType = 'PBS'; 
} 
my $schedulerDir = qx|which qstat|;
$schedulerDir or $schedulerDir = "";
chomp $schedulerDir;
$schedulerDir =~ s|/qstat||;
#------------------------------------------------------------------------
sub isShellCommand {
    my ($shellCommand) = @_;
    return !system($bashPath, "-c", "which $shellCommand &> /dev/null");
}
#------------------------------------------------------------------------
# parse the path to the q program and environment targets
#------------------------------------------------------------------------
print ".";
$script = "$qDir/q";
my $libDir = "$qDir/lib";
my $modulesDir = "$qDir/modules";
my $utilitiesDir = "$qDir/utilities";
my $utilitiesPath = "";
foreach my $dir(<$utilitiesDir/*>){
    -d $dir or next;
    $utilitiesPath .= "$dir:";
}
my $envScript = "$qDir/remote/def_env.pl";
#------------------------------------------------------------------------
# populate list of known system commands
#------------------------------------------------------------------------
my %execs;
foreach my $path(split(":", "$utilitiesPath$ENV{PATH}")){
    print ".";
    foreach my $file(<$path/*>){
        -d $file and next;
        -x $file or next;
        $file =~ s|^$path/||;
        $execs{$file}++;
    }
}
my $execs = join(" ", keys %execs);
#------------------------------------------------------------------------
# print the q program target script
#------------------------------------------------------------------------
print "\ngenerating q program target\n";
open my $outH, ">", $script or die "could not open $script for writing: $!\n";
print $outH
'#!'.$perlPath.'
use strict;
use warnings;
$ENV{PATH} = "'.$utilitiesPath.'$ENV{PATH}";
our $version = '."'$version'".';
our $perlPath = '."'$perlPath'".';
our $bashPath = '."'$bashPath'".';
our $timePath = '."'$timePath'".';
our $timeVersion = '."'$timeVersion'".';
our $memoryCorrection = '."'$memoryCorrection'".';
our $memoryMessage = '."'$memoryMessage'".';
our $qDir = '."'$qDir'".';
our $libDir = '."'$libDir'".';
our $modulesDir = '."'$modulesDir'".';
our $qType = '."'$qType'".';
our $schedulerDir = '."'$schedulerDir'".';
$ENV{Q_Q_TYPE} = '."'$qType'".';
$ENV{Q_UTIL_DIR} = '."'$utilitiesDir'".';
$ENV{Q_MOD_DIR} = '."'$modulesDir'".';
our %shellCommands = map { $_ => 1 } qw('.$execs.');
$| = 1;
require "$qDir/lib/main.pl";
qMain();
';
close $outH;
#------------------------------------------------------------------------
# make the q program target script executable
#------------------------------------------------------------------------
qx|chmod ugo+x $script|;
#------------------------------------------------------------------------
# collect default environment variables for use by q remote server mode
#------------------------------------------------------------------------
print "storing default environment information\n";
open $outH, ">", $envScript or die "could not open $envScript for writing: $!\n";
print $outH "use strict;\n";
print $outH "use warnings;\n";
my %ignore = map { $_ => 1 } qw(HOME HOSTNAME LOGNAME MAIL USER USERNAME);
foreach my $key(keys %ENV){
    $ignore{$key} and next;
    print $outH "\$ENV{'$key'} = '$ENV{$key}';\n";    
}
print $outH "1;\n";
close $outH;
#------------------------------------------------------------------------
# copy the version file for use by q remote
#------------------------------------------------------------------------
system("cp $versionFile $qDir/remote");
#------------------------------------------------------------------------
print "done\n";
print "created q program target:\n  $script\n";
#========================================================================

