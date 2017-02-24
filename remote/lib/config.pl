#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($scriptDir %options $isServer $debug);
our %config = ();
our $configFile = $isServer ? "$scriptDir/server.conf" : "$scriptDir/$options{config}";
my @requiredOptions = $isServer ? qw(userKey qCommand) : qw(host); 
my %optionDefaults = (localPort    => 8000,  #keyFile default deterimined by ssh
                      sshCommand   => 'ssh',
                      scpCommand   => 'scp',
                      agentCommand => 'ssh-add',
                      port         => 22,  # directory default determined at login      
                      user         => $ENV{USER} ? $ENV{USER} : $ENV{USERNAME},
                      qCommand     => 'q');
############################################################

############################################################
# load user configuration
#-----------------------------------------------------------
sub loadConfig{  # read q_remote.conf
    -e $configFile or die "could not find configuration file $configFile\n";
    my $config = slurpFile($configFile);
    $config = stripComments(\$config); 
    loadConfigLines(\$config, \%config);        
    checkConfig();     
}
sub loadConfigLines {  # push config lines into hashes
    my ($lines, $hash) = @_;
    open my $inH, "<", $lines;
    while(my $line = <$inH>){
        chomp $line;
        $line or next;
        $line =~ s/^\s+//g;
        $line =~ s/\s+$//g;
        my ($param, @value) = split(/\s+/, $line);
        my $value = join(" ", @value);
        $$hash{$param} = $value;
    }
    close $inH;
}
sub checkConfig {  # fill in system configuration defaults, check for option errors
    foreach my $option(keys %optionDefaults){ $config{$option} ||= $optionDefaults{$option} }
    foreach my $option(@requiredOptions){ $config{$option} or  die "config error: '$option' is required\n" }
    if($isServer){
        $config{user} = $ENV{$config{userKey}};
        $config{user} or die "config error: \$ENV\{$config{userKey}\} has no value\n";
    } else {
        $config{localPort} =~ m|\D| and die "config error: invalid localPort specified'\n";
        $config{port} =~ m|\D| and die "config error: invalid port specified'\n";
        $config{portOption} = $config{port} != 22 ? "-p $config{port}" : "";   
        $config{batchOption} = "";
        lc($config{sshCommand}) =~ m|plink| and $config{batchOption} = '-batch';    
        $config{sshCommand} =~ m|ssh$| and $config{batchOption} = '-o "BatchMode yes"';
        $config{keyFileOption} = $config{keyFile} ? "-i $config{keyFile}" : "";
    }
}
############################################################

1;

