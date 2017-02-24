#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw(%config $params %project $isServer $debug);
our $isSSH = undef;
our $sshError = undef;
my %chmod = map { $_ => 1 } qw(submit status);
############################################################

############################################################
# interface with the execution host
#-----------------------------------------------------------
sub checkSSH {  # check logon on first encounter
    $sshError = "";
    $isServer and $isSSH = 1;
    $isSSH and return;  # is a known SSH-capable configuration
    my $OK = "user home is: ";
    my $ssh = getSSH("echo $OK\$HOME");
    $ssh =~ m|$OK(.+)| and $isSSH = 1;
}
sub getSSH {  # execute a command on remote host and return result en bloc
    my ($commandString) = @_;
    if($isServer){
        return qx|$commandString|;
    } else {
        my $remote = getRemote($commandString);
        return qx|$remote|;  # open "-|" and ssh/plink -t option allows line-by-line retrieval
                             # and HTTP::Daemon can use callback to submit lines one at a time
                             # HOWEVER, browsers invariably buffer on receipt!
    }
}
sub getRemote {  # parse the ssh string
    my ($commandString, $useTTY) = @_;
    my $tty = $useTTY ? '-t' : '';
    # NOTE: getRemote does NOT automatically single- or double-quote commandString
    # however, quoting IS required if command string contains concatenated commands
    # of form "command1 ; command2".  In this case, if commandString is not quoted
    # by the calling sub, command1 is passed to ssh, but command2 is executed on  
    # daemon host, not on the remote server as is probably intended
    return "$config{sshCommand} $config{batchOption} $tty ".
           "$config{portOption} $config{keyFileOption} ".
           "$config{user}\@$config{host} $commandString";
}
sub scpUpload {  # copy a file to host
    my ($localFile, $remoteFile) = @_;
    my $remote = "$config{scpCommand} $config{batchOption} ".
                 "$config{portOption} $config{keyFileOption} ".
                 "$localFile $config{user}\@$config{host}:$remoteFile";
    qx|$remote|;
}
sub scpDownload {  # copy a file from host
    my ($localFile, $remoteFile) = @_;
    my $remote = "$config{scpCommand} $config{batchOption} ".
                 "$config{portOption} $config{keyFileOption} ".
                 "$config{user}\@$config{host}:$remoteFile $localFile";
    qx|$remote|;
}
sub getQH {  # execute q on remote host and return handle to result stream
    my (@qRefs) = @_;
    my $masterDir = "$project{directory}/masters/$$params{masterClass}";
    my $master = "$masterDir/$$params{masterName}";
    my $qDataDir = "$masterDir/.$$params{masterName}.data";
    my @qs;
    foreach my $qRef(@qRefs){
        my ($qCommand, @options) = @$qRef;
        $$params{isDryRun} and push @options, '--dry-run'; 
        push @options, "--_q_remote_"; 
        $isServer and push @options, "--_server_mode_"; 
        my $options = join(" ", @options);
        my $q = "$config{qCommand} $qCommand $options $master 2>&1";  
        push @qs, $q;
        $isServer and $chmod{$qCommand} and push @qs, "chmod -R ug+rw $qDataDir";
    }
    my $q = join(" ; ", @qs);  # allow composite command sequences (e.g. report, script, environment)
    $q = $isServer ? $q : getRemote('"'.$q.'"', 1);
    $q or return; 
    open my $qH, "-|", $q or ($debug .= "$!<br>\n" and return);
    return $qH;  # calling sub must close $qH
}
############################################################

1;

