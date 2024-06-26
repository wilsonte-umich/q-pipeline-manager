
# This file contains functions useful to Perl worker scripts.
# Initialize them in your worker script by adding line:
#
# require "$ENV{Q_UTIL_DIR}/utilities.pl";

#####################################################################
# functions that interpret array job $TASK_ID
#--------------------------------------------------------------------
# usage:
#     my $TASK_ID = getTaskID();
#     my $taskObject = getTaskObject("LIST_NAME");  # LIST_NAME refers to a space-delimited environment variable
#--------------------------------------------------------------------
sub getTaskID {
    my $TASK_ID = $ENV{SGE_TASK_ID} ? $ENV{SGE_TASK_ID} : ($ENV{PBS_ARRAYID} ? $ENV{PBS_ARRAYID} : $ENV{SLURM_ARRAY_TASK_ID});
    defined $TASK_ID or die "getTaskID error: neither SGE_TASK_ID nor PBS_ARRAYID nor SLURM_ARRAY_TASK_ID was set\n".
                            "most likely, job scheduler option -t was not properly configured\n";
    return $TASK_ID;
}
#--------------------------------------------------------------------
sub getTaskObject {
    my ($LIST_NAME) = @_;  # e.g. getTaskObject("SAMPLES")
    defined $LIST_NAME or die "getTaskObject error: missing LIST_NAME\n";
    defined $ENV{$LIST_NAME} or die "getTaskObject error: environment variable $LIST_NAME not set\n";
    my $TASK_ID = getTaskID();    
    my $i = $TASK_ID - 1;
    my @list = split(/\s+/, $ENV{$LIST_NAME});
    my $taskObject = $list[$i];
    defined $taskObject or die "getTaskObject error: environment variable $LIST_NAME did not have element $i\n";
    return $taskObject;
}
#####################################################################


#####################################################################
# pipeline file handling functions
#--------------------------------------------------------------------
# usage:
#     waitForFile($FILE_NAME[,$DIE_ON_FAIL,$TIME_OUT]);  # default time out = 60 seconds; $DIE_ON_FAIL = false
#--------------------------------------------------------------------
sub waitForFile {  # wait for expected file from previous pipeline step to appear on file system
    my ($FILE_NAME, $DIE_ON_FAIL, $TIME_OUT) = @_;  # e.g. waitForFile("/path/to/my/file")
    $FILE_NAME or die "waitForFile error: missing FILE_NAME\n";
    $TIME_OUT or $TIME_OUT = 60;
    my $elapsed = 0;
    do { -e $FILE_NAME and return; sleep 1; $elapsed++ } until ($elapsed > $TIME_OUT);
    $DIE_ON_FAIL and die "waitForFile error: $FILE_NAME not found in $TIME_OUT seconds\n";
}
#####################################################################

1;

