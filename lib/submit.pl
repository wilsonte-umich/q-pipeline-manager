use strict;
use warnings;
use Cwd(qw(getcwd abs_path));

#========================================================================
# 'submit.pl' is the q instructions file interpreter and qsub caller
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($qDir $qType $schedulerDir %options %optionInfo %commandOptions   
            %protectFileGlobs %backupDirs @instrsTree $memoryMessage $bashPath $timePath
            $masterFile $masterDir $masterFileName $statusFile $logDir $modulesDir);
my ($runDir, $ovrideJobName, $ovrideRunDir, $isJobNameTag, $isRunDirTag,
    %variables, %lastJobInThread, %threadPred, $dependencies, %directives,
    %jobIDs, @statusInfo, %mergeThreads, %returnVariables, %jobInfos,
    %masterThreadNumbers, $jobsAdded, $needAutoProtect, $needAutoBackup,
    $shebang, $isBash, $qsubMode, $qsubScriptFile, %currentInvokeVars,
    $currentInstrsDir, %fileContents, %fileModuleCache, $embeddedTreeFile);
our ($jobName, $currentJob, $array, %jobPred, %jobSucc, $getOutputFiles,  
     $qsubOptions, $directives, $qInUse, %currentInstrsFile, $currentScriptFile, %virtualFiles);
our $threadLevel = 0; #the root master level
our $threadCounter = 0; #implied first thread of master instructions file
our %currentThread = ($threadLevel => $threadCounter);
my $qsubInvokeTarget = '__qsub_invoke_target__';
my $qTarget = abs_path($0);
my $fileModulePath = "$modulesDir/file";
my %cacheSlaves = ("$fileModulePath/exists.q" => ['FILE', ['EXISTS', 'N_FILES']], 
                   #"$fileModulePath/require.q" => ['FILE', ['N_FILES']],  
                   "$fileModulePath/gzipped.q" => ['FILE', ['GZIPPED']],  
                   "$fileModulePath/check_prefix.q" => ['PREFIX', []], 
                   "$fileModulePath/create.q" => ['DIR', []]);  
my %conditionals = map { $_ => 1 } qw(exitIf exitUnless dieIf dieUnless);   
#========================================================================

#========================================================================
# main execution block
#------------------------------------------------------------------------
sub qSubmit {      
    checkScheduler();
    checkLock();
    checkDeleteExtend();
    checkSyntax('submit');
    followInstructions();
    setAutoUnprotect();
    provideFeedback();  
}
#------------------------------------------------------------------------
sub checkScheduler {  # make sure there will be a way to run requested jobs
    $qInUse = $options{'execute'} ? 'local' : $qType;
    $qInUse or die "job submission requires either Sun Grid Engine (SGE), Torque/PBS, or option --execute\n";
}
sub checkSyntax {  # check syntax in dry-run before submitting any jobs
    my ($command) = @_;
    $options{'dry-run'} and return;  # only run syntax check first if this is an execution run 
    print "checking pipeline syntax\n";
    my $options = '--dry-run --_suppress-echo_';
    foreach my $option(keys %options){  # pass current options to syntax check run 
        defined $options{$option} or next;
        my $value = ${$optionInfo{$option}}[1] ? $options{$option} : "";
        $options .= " --$option $value" 
    } 
    my $syntaxCheck = "perl $qTarget $command $options $masterFile"; # syntax check executed by recalling q
    my $syntaxError = system($syntaxCheck); # returns zero if syntax check run of q does not die
    $syntaxError and die "unrecoverable error, no jobs queued\n";  #main run dies on syntax error
}
sub followInstructions {  # execute jobs in the required order 
    $options{'_suppress-echo_'} or print "job_name                    array  job_ID   job_#  depends on job_#\n";
    executeInstructions( getEnvFiles('q') );  # environment files always acted on first
    initializeAutoThread('pre');
    my $needAutoUnprotect = needAutoUnprotect();
    addAutoJob($needAutoUnprotect, 'unprotect', '500', '1:00:00', '-w u'); # remove write protection *before* running user jobs 
    initializeAutoThread('main');
    executeInstructions($masterFile); 
    initializeAutoThread('post');
    addAutoJob($needAutoProtect, 'protect', '500', '1:00:00'); # add write protection *after* running user jobs
    addAutoJob($needAutoBackup, 'backup', '500', '10:00:00');  # add output archiving *after* running user jobs
    addAutoJob($options{publish}, 'publish', '500', '1:00:00', getPassedAutoOptions('publish')); # publish pipeline *after* running user jobs  
}
sub provideFeedback {  # exit feedback
    if ($options{'dry-run'}){
        print "no syntax errors detected\n";
    } elsif($jobsAdded) {  
        generateStatusFile(); # generate disk copy of queued jobs
        print "all jobs queued\n";  
    } else {
        print "no jobs to queue\n";  
    } 
}
#========================================================================

#========================================================================
# find instructions files
#------------------------------------------------------------------------
sub getEnvFiles { # environment files, in order of increasing precedence
    my ($suffix) = @_;
    my $fileRoot = "environment.$suffix";
    return ("$qDir/$fileRoot", 
            "$ENV{HOME}/.q/$fileRoot",
            "$masterDir/../$fileRoot",
            "$masterDir/$fileRoot")
}
sub findInstrsFile { # attempt to find a slave instructions file
    my ($instrsFile) = @_;
    $virtualFiles{$instrsFile} and return \$virtualFiles{$instrsFile}; # virtual files declared in master (or slave) take precedence over disk files
    my @dirs =  ($masterDir,         # look first in same directory as the master instructions file
                 $currentInstrsDir,  # then look in the same directory as the invoking instructions file
                 $modulesDir,        # then look in the q modules directory
                 '',                 # then attempt to interpret as relative or absolute file path
                 );       
    my @searchedDirs;
    foreach my $i(0..$#dirs){
        defined $dirs[$i] or next; 
        $dirs[$i] and push @searchedDirs, $dirs[$i];       
        my $file = $dirs[$i] ? "$dirs[$i]/$instrsFile" : $instrsFile;  
        if($fileContents{$file} or -e $file){
            $file =~ m|^\/| or $file = abs_path($file); 
            return $file;
        }
    }
    my $searchedDirs = join("\n  ", @searchedDirs);
    die "could not find instructions file $instrsFile in:\n  $searchedDirs\n";
}
#========================================================================

#========================================================================
# parse instructions files
#------------------------------------------------------------------------
sub executeInstructions { # parse and execute the instructions in instrsFile(s)
    my (@instrsFiles) = @_;
    foreach my $instrsFile(@instrsFiles){
        my $instrs;
        if(ref($instrsFile) eq 'SCALAR'){
            push @instrsTree, [$threadLevel, $embeddedTreeFile, "embedded"];
            $currentInstrsFile{$threadLevel} = $currentInstrsFile{$threadLevel-1}; 
            $instrs = slurpFile($instrsFile); 
        } elsif($fileContents{$instrsFile}){
            push @instrsTree, [$threadLevel, $instrsFile];
            $currentInstrsFile{$threadLevel} = $instrsFile;
            $instrs = $fileContents{$instrsFile};
            ($currentInstrsDir) = $instrsFile =~ m|(.+)/.+|;
        } elsif(-e $instrsFile){
            push @instrsTree, [$threadLevel, $instrsFile];
            $currentInstrsFile{$threadLevel} = $instrsFile;
            $instrs = slurpFile($instrsFile); 
            $fileContents{$instrsFile} = $instrs;  # maintain a cache of encountered instructions files
            ($currentInstrsDir) = $instrsFile =~ m|(.+)/.+|;   
        } else {
            next;        
        }
        my %threadNumbers;  # store conversion of user thread names to thread numbers on a per-file basis  
        $instrsFile eq $masterFile and %threadNumbers = %masterThreadNumbers;   
        while($instrs =~ s|<file\s+name="{0,1}(.+?)"{0,1}>(.*?)</file>||s){ $virtualFiles{$1} = $2 } # capture and remove any virtual files
        open my $instrsFileH, "<", \$instrs;
        while (my $line = <$instrsFileH>){
            defined processInstruction($line, $instrsFile, $instrsFileH, \%threadNumbers) or last;
        } 
        close $instrsFileH;
        $instrsFile eq $masterFile and %masterThreadNumbers = %threadNumbers;  # remember master thread numbers for auto jobs
    }
}
sub processInstruction { # main interpreter switch
    my ($line, $instrsFile, $instrsFileH, $threadNumbers) = @_;
    my $parsedLine = parseInstructionLine($line, $instrsFileH);
    $parsedLine or return 0;
    my $genericLine = $parsedLine;
    $parsedLine = replaceVariables($parsedLine);
    my @line = split(/\s/, $parsedLine);
    my $instruction = $line[0];
    my $currentThread = getCurrentThread($threadLevel);
    if ($instruction =~ m/^\$(.+)/){ # variable assignment lines
        assignVariable($parsedLine, $1, @line[1..$#line]);
    } elsif ($instruction eq 'exit'){ # exit control from q files
        return;
    } elsif ($instruction eq 'exitUnless'){
        $getOutputFiles and return 1;  # ignore flow control when only collecting file targets
        checkExitCondition($parsedLine, $instruction, @line[1..$#line]) or return;
    } elsif ($instruction eq 'exitIf'){
        $getOutputFiles and return 1;
        checkExitCondition($parsedLine, $instruction, @line[1..$#line]) and return;
    } elsif ($instruction eq 'dieUnless'){
        $getOutputFiles and return 1;
        checkExitCondition($parsedLine, $instruction, @line[1..$#line]) or die "died at line:\n$parsedLine\n";
    } elsif ($instruction eq 'dieIf'){
        $getOutputFiles and return 1;
        checkExitCondition($parsedLine, $instruction, @line[1..$#line]) and die "died at line:\n$parsedLine\n"; 
    } elsif ($instruction eq 'thread'){ # parallelization commands
        $threadNumbers or die "'#q thread' not allowed as target script directive\n";
        $line[1] or die "missing thread name in $genericLine\n";
        $instrsFile eq $masterFile and push @line, 'q_auto_pre';  #all declared master threads depend on thread q_auto_pre
        initializeThread($instrsFile, $threadNumbers, @line);
        $options{'_suppress-echo_'} or ($options{'verbose'} and print "@line\n");        
    } elsif ($instruction eq 'invoke'){ # slave calls
        my $firstOfSeries = 1;
        processInvoke($currentThread, \$firstOfSeries, @line);
    } elsif ($instruction eq 'preserve'){ #slave return values
        preserveVariables($parsedLine, $instrsFile, @line);
    } elsif ($instruction eq 'qsub'){ # queue target jobs
        if($qsubMode == 4) {  # handle qsub matrix invoke on behalf of user
            my $qSF = $qsubScriptFile;
            $qSF =~ s/\$/\\\$/g;
            $line =~ s/qsub\s+$qSF/invoke $qsubInvokeTarget/;
            processInstruction($line, $instrsFile, $instrsFileH, $threadNumbers, 1);
        } else {  # single qsub call, handle immediately
            queueScript($currentThread, $genericLine, @line[1..$#line]);
            if($currentScriptFile){
                push @instrsTree, [$threadLevel + 1, $currentScriptFile];
            } else {
                push @instrsTree, [$threadLevel + 1, $line[1], "embedded"];
            }
        }  
    } elsif ($instruction eq 'protect'){
        addAutoFileGlob($genericLine, \%protectFileGlobs, @line[1..$#line]); 
        $needAutoProtect = 1;
    } elsif ($instruction eq 'backup'){
        addAutoDirectory($genericLine, \%backupDirs, @line[1..$#line]); 
        $needAutoBackup = 1;
    } elsif ($instruction eq 'backupDir'){
        $options{'_backup-dir_'} and print "CAUTION: only the last encountered 'backupDir' instruction is used!!\n";
        $options{'_backup-dir_'} = $line[1];
    } elsif ($instruction eq 'publishDir'){
        $options{'out-dir'} = $line[1];
    } elsif ($instruction eq 'publishTitle'){
        $options{'title'} = join(" ", @line[1..$#line]);
    } elsif ($instruction eq 'publishIntro'){
        $options{'intro-file'} = abs_path($line[1]); 
    } elsif ($instruction eq 'publishMask'){
        my $masks = join(",", @line[1..$#line]);
        $options{'mask'} = $options{'mask'} ? "$options{'mask'},$masks" : $masks;
    } elsif(isShellCommand($instruction)){ # commands to be run in shell at submission time
        $getOutputFiles and return 1;
        runImmediately(undef, @line);
    } else { # unrecognized by q or shell
        $options{'_suppress-echo_'} or print "\n";
        die "'$instruction' is not a recognized parameter or command in line:\n  $parsedLine\n";
    }
    return 1;  
}
my $inBlock;
sub parseInstructionLine { # deconstruct the complexity of the typed instruction line
	no warnings 'recursion';
    my ($line, $instrsFileH) = @_;
    $line or return "";
    chomp $line;
    $line =~ s/\r//g;
    $line =~ m/^\s*([^#]*)/; # strip leading white space and ignore comments
    $1 or return ""; # ignore lines with no content
    $line = $1;    
    if($line =~ m/^>+\s*$/){ $inBlock = 1; $line = "" } # obey line block demarcation
    if($line =~ m/(.*)\\\s*$/){ # recursively append continuation lines, bash style
        $line = $1;
        $line = appendNextLine($line, $instrsFileH);
    } elsif ($inBlock){ # >>> <<< style
        $line = appendNextLine($line, $instrsFileH);
    }    
    $line =~ s/\t/ /g; # all lines collapsed to single-space delimited word series
    $line =~ tr/ / /s;      
	$line =~ s/\s+$//; # strip trailing white space
	$line =~ s/^\s+//; # strip leading white space again (added by >>> <<< blocks)
    return $line; 
}
sub appendNextLine { # handle multi-line q instructions
	no warnings 'recursion';
    my ($line, $instrsFileH) = @_;
    $instrsFileH or die "continuation lines not allowed in target script directives\n";
    my $parsedNextLine;
    until($parsedNextLine){
        my $nextLine = <$instrsFileH>;
        $nextLine or last;
        if($nextLine =~ m/^<+\s*$/){ $inBlock = undef; last}
        $parsedNextLine = parseInstructionLine($nextLine, $instrsFileH);
    }
    $parsedNextLine and $line = "$line $parsedNextLine";   
    return $line;
}
#========================================================================

#========================================================================
# instruction variable handling
#------------------------------------------------------------------------
sub replaceVariables { # replace variable calls with most current assignment values
    my ($line, $isCircuit) = @_;
    $line =~ m/^(\S+) (.+)$/ or return $line; # nothing to replace
    my ($instruction, $arguments) = ($1, $2);
    $arguments =~ s/\\\$/~masked~escaped~variable~/g; # mask escaped variables
    $arguments =~ s/\$NULL/~masked~variable~NULL/g;   # mask nulls
    if ("\U$instruction" eq 'INVOKE'){                # mask invocation slave variables
        maskInvocationMatrix(\$arguments);
    } elsif ("\U$instruction" eq 'QSUB'){             # mask invocation slave variables in qsub if needed  
        setQsubMode($line, split(/\s/, $arguments));
        $qsubMode == 4 and maskInvocationMatrix(\$arguments);
    } elsif ("\U$instruction" eq 'PRESERVE'){         # mask slave return variables
        $arguments =~ s/^\$\$/~masked~reference~/g;    
        $arguments =~ s/^\$/~masked~variable~/g;
    }
    my @varChunks = split(/\$/, $arguments);
    my $leader = shift @varChunks;
    scalar(@varChunks) or return $line; # nothing to replace
    foreach my $i(0..$#varChunks){
        $varChunks[$i] =~ m/^(\w+)(.*)$/ or next;
        my ($var, $rest) = ($1, $2);
        unless(defined $variables{$var}){
            if($isCircuit){
                $variables{$var} = '__circuit_filler_value__';  
            } elsif($i == 0 and $conditionals{$instruction}){  # allow optional variables in exit conditions
                $variables{$var} = '$NULL';
            } elsif($i == 0 and $instruction =~ m|^\$| and $line =~ m/ \? .+ \: /){  # allow optional variables as ternary operator test case
                $variables{$var} = '$NULL';
            } else {
                die "variable \$$var used before being assigned a value in line:\n$line"
            }
        }
        defined $rest or $rest = '';
        $varChunks[$i] = "$variables{$var}$rest";
    }
    $arguments = join("", $leader, @varChunks);
    $arguments =~ s/~masked~escaped~variable~/\\\$/g; # unmask
    $arguments =~ s/~masked~reference~/\$\$/g;    
    $arguments =~ s/~masked~variable~/\$/g;
    return "$instruction $arguments";
}
sub maskInvocationMatrix {
    my ($arguments) = @_;
    $$arguments =~ s/^(\S+ )\$/$1~masked~variable~/;
    $$arguments =~ s/( \+ )\$/$1~masked~variable~/g;
    $$arguments =~ s/( \* )\$/$1~masked~variable~/g;
}
sub assignVariable { # handle assignment of values to variables
    my ($parsedLine, $varName, @valueArray) = @_;
    $varName eq 'NULL' and die "cannot redefine variable \$NULL\n";
    if($parsedLine =~ m/^\$$varName.* \?.+\: /){ # allow variable assignment via ternary operator
        @valueArray = applyTernaryOperator($parsedLine, @valueArray);
    } elsif ($parsedLine =~ m/^\$$varName run/ or # allow variable assignment from shell command result
             $parsedLine =~ m/^preserve \$$varName run/){
        @valueArray = (runImmediately(1, @valueArray[1..$#valueArray]));
    } 
    fillVariable($varName, join(" ", @valueArray));
    $options{'_suppress-echo_'} or ($options{'verbose'} and print "\$$varName $variables{$varName}\n");   
}
sub fillVariable {
    my ($varName, $varValue) = @_;
    $variables{$varName} = $varValue;
    "\U$variables{$varName}" eq 'FALSE' and $variables{$varName} = '$NULL';      
    $ENV{$varName} = $variables{$varName}; #mimic all varValues to ENV for passing to target script via -V
    $ENV{$varName} =~ s/\$NULL//g;    # remove null placeholder
    $ENV{$varName} =~ s/\\(\\*)/$1/g; # remove q's variable escape character   
}
sub clearVariable {
    my ($varName) = @_;
    delete $variables{$varName};
    delete $ENV{$varName};
}
sub resetVariable {
    my ($varName, $resetValues) = @_;
    $variables{$varName} = $$resetValues{$varName};
    $ENV{$varName} = $variables{$varName};
}
sub applyTernaryOperator { # interpret ternary operator conditional: $VAR test_case ? true_value : false_value
    my ($parsedLine, @ternary) = @_;
    my $die = "incorrect use of ternary operator in line $parsedLine\nusage = \$VAR test_case ? true_value : false_value\n";
    my ($i1, $i2); # index of ? and : operators, respectively
    foreach my $i (0..$#ternary){
        $ternary[$i] or next;
        if($ternary[$i] eq '?'){
            defined $i1 and die $die; # must only be one instance of each operator
            $i1 = $i;
        } elsif ($ternary[$i] eq ':'){
            defined $i2 and die $die;
            $i2 = $i;
        }
    }
    ($i1 and $i2) or die $die; 
    $i2 > $i1 + 1 or die $die; # ensure provision of values between operators
    my $testCase = join(" ", @ternary[0..$i1-1]); # assemble parts to allow case and values to be multi-word assignments
    my $trueValue = join(" ", @ternary[$i1+1..$i2-1]);
    my $falseValue = join(" ", @ternary[$i2+1..$#ternary]);
    if(defined $testCase and !($testCase eq '' or $testCase eq '$NULL')){
        return getTernaryResult($trueValue);    
    } else {
        return getTernaryResult($falseValue);
    }
}
sub getTernaryResult {
    my ($value) = @_;
    if(defined $value){
        $value =~ m/^run (.+)/ and $value = runImmediately(1, split(" ", $1));
        return $value;
    } else {
        return '$NULL';
    }
} 
sub checkExitCondition {
    my ($parsedLine, $instruction, @valueArray) = @_;
    assignVariable("\$$parsedLine", $instruction, @valueArray);
    return ($variables{$instruction} and $variables{$instruction} ne '$NULL');
}    
sub preserveVariables {
    my ($parsedLine, $instrsFile, @line) = @_;
    ($line[1] =~ m/^\$(.+)/ or $line[1] =~ m/^(all)\s*$/) or 
        die "could not extract variable name from preserve statement:\n$parsedLine\n";
    my $varName = $1;        
    if($varName eq 'all'){
        foreach my $varName(keys %variables){ 
            $currentInvokeVars{$varName} and die "'preserve' cannot use master variable '\$$varName'; it is already in use as an invoke variable\n";
            push @{$returnVariables{$instrsFile}}, [$varName];
        }
    } else {
        my $isReference;
        if($varName =~ m/^\$(.+)/){  # master variable is in format '$$varName'
            $varName = $variables{$1};  # preserve into master variable whose name is given by the value of $varName        
            $varName or die "could not recover master variable name from line:\n $parsedLine\n";
            $varName =~ m|\W| and die "master variable reference contained character(s) that cannot be used in a variable name:\n $parsedLine\n";
            $isReference = 1;
        }
        $currentInvokeVars{$varName} and die "'preserve' cannot use master variable '\$$varName'; it is already in use as an invoke variable\n";
        assignVariable($parsedLine, $varName, @line[2..$#line]);
        push @{$returnVariables{$instrsFile}}, [$varName, $isReference];
    }
}
#========================================================================

#========================================================================
# thread management
#------------------------------------------------------------------------
sub getCurrentThread { # determine the current thread for a level; may come from a master if a slave declares no threads
    my ($threadLevel) = @_;
    my $i = $threadLevel; # loop will always terminate at threadLevel=0, if not before
    until (defined $currentThread{$i}) { $i-- } 
    return $currentThread{$i};
}
sub initializeThread { # thread calls break the dependency chain to all except declared merge threads
    my ($instrsFile, $threadNumbers, $threadKeyword, $threadName, @mergeThreads) = @_;
    $$threadNumbers{$threadName} and die "thread $threadName redefined in $instrsFile\n"; #threads name must be unique within a file 
    $threadCounter++; # threads are simply numbered in order for tracking
    my $currentThread = $threadCounter;
    $currentThread{$threadLevel} = $currentThread; 
    $$threadNumbers{$threadName} = $currentThread; 
    foreach my $i(0..$#mergeThreads){  # merge threads must be previously named within the same file
        my $mergeThread = $$threadNumbers{$mergeThreads[$i]};  # replace user-supplied thread names with internal thread numbers
        defined $mergeThread or die "thread $mergeThreads[$i] was called as a dependency before it was defined in $instrsFile\n";
        my @ljit = getLastJobsInThread($mergeThread);
        scalar(@ljit) and push @{$threadPred{$currentThread}}, @ljit; # store thread merge dependencies for first job in this thread
        $mergeThreads[$i] = $mergeThread; # convert threads names to number for saving in %mergeThreads
    }
    $mergeThreads{$currentThread} = \@mergeThreads; # remember merge threads for later deletion of master thread dependencies on them
}
sub getLastJobsInThread { # the job or jobs queued by the most recent line in a thread; can be multiple when slaves execute many threads
    my ($thread, $invoke) = @_; # $invoke is a flag to return the values being held by an invoking thread
    defined $lastJobInThread{$thread} or return;  
    my $jobs = $invoke ? $lastJobInThread{$thread}{invoke} : $lastJobInThread{$thread};    
    my @ljit;
    getLastJobsInThread_($jobs, \@ljit);
    return @ljit;
}
sub getLastJobsInThread_ { # every thread/slave chain ultimately ends in a single queued job; recursively find and collect these
    my ($jobs, $ljit) = @_;
    if(ref($jobs) eq "HASH"){ 
        foreach my $key(keys %$jobs){ getLastJobsInThread_($$jobs{$key}, $ljit) }
    } else {
        push @$ljit, $jobs; 
    }
}
#========================================================================

#========================================================================
# slave invocation
#------------------------------------------------------------------------
sub processInvoke { # handling invoke calls to slaves
    my ($currentThread, $firstOfSeries, @line) = @_;
    my $slaveFile = $line[1];
    
    $slaveFile =~ s/\\(\\*)/$1/g; # remove q's variable escape character
    
    my %tmp;
    $lastJobInThread{$currentThread} and %tmp = %{$lastJobInThread{$currentThread}};
    delete $lastJobInThread{$currentThread}; # clear current thread job chain in preparation for slave job additions
    %{$lastJobInThread{$currentThread}{invoke}} = %tmp; # hold current values for determining invoked thread dependence
    my %tmpInvokeVars = %currentInvokeVars;
    %currentInvokeVars = ();
    if($line[2]){ # variable assignment(s) present; is matrix invoke
        my %priorVars = %variables;  # remember variables as they existed prior to invoke in case a master variable is used as an invoke variable
        extractSlaveMatrix(\@line, \my@slaveVars);  
        foreach my $slaveVar(@slaveVars){
            my ($varName) = @$slaveVar;
            $currentInvokeVars{$varName} = 1;  # used to ensure that a slave doesn't use an invoke variable as a preserve master variable
        }  
        invokeSlaveMatrix($firstOfSeries, $slaveFile, \@slaveVars, 0);
        foreach my $slaveVar(@slaveVars){  # reset invoke variables to whatever status they had in master prior to invoke
            my ($varName) = @$slaveVar;
            exists $priorVars{$varName} ? resetVariable($varName, \%priorVars) : clearVariable($varName);
        }       
    } else { # simple invoke using only prior out-of-statement variable assignments
        invokeSlave($firstOfSeries, $slaveFile);
    }     
    %currentInvokeVars = %tmpInvokeVars;
    delete $lastJobInThread{$currentThread}{invoke}; # clear held values
    if(scalar(keys %{$lastJobInThread{$currentThread}})){ # invocation added new job dependencies
        updateMasterDependencies($threadLevel, $currentThread); # always update masters when slaves change
    } else { # was a dry invoke, e.g. only changed parameters
        %{$lastJobInThread{$currentThread}} = %tmp; # preserve prior dependencies when slaves do not queue any jobs
    } 
}
sub extractSlaveMatrix { # handles matrix invoke pragma: invoke my.slave $VAR1 $VAR1S : $VAR2 $VAR2S
    my ($line, $slaveVars) = @_;
    my $i = 2; # enter loop on first slaveVar
    while($$line[$i]){
        $$line[$i] =~ m/^\$(.+)/ or die "error in matrix invoke: $$line[$i] is not a variable name\n";
        my $slaveVar = $1;
        $i++; # move to first assignment value        
        my @masterValues;
        while (defined $$line[$i] and !($$line[$i] eq '+' or $$line[$i] eq '*')){ # collect list of master assignments to slaveVar
            push @masterValues, $$line[$i]; 
            $i++;
        } 
        push @$slaveVars, [$slaveVar, $$line[$i], @masterValues]; # collect slaveVars and assigned master values
        $i++; # move to next slaveVar, if any    
    }
}
sub invokeSlaveMatrix { # commit all possible combinations of master assignments to slave variables with repeated invocation of slave
    my ($firstOfSeries, $slaveFile, $slaveVars, $i) = @_;
    if($$slaveVars[$i]){
        my ($slaveVar, $matrixType, @masterValues) = @{$$slaveVars[$i]};
        foreach my $j(0..$#masterValues){ 
            fillVariable($slaveVar, $masterValues[$j]); # assign master values to slaveVar one at a time
            my ($k, $mType, $slVar, @masterVals) = (0, $matrixType);
            while($mType and $mType eq '+'){ # apply linear matrix by simultaneously committing jth assigned value of all variables
                $k++;
                $$slaveVars[$i + $k] or die "'\$$slVar' followed by linear matrix operator '+' but subsequent variable is missing\n";
                ($slVar, $mType, @masterVals) = @{$$slaveVars[$i + $k]}; # values of next slaveVar in linear matrix
                defined $masterVals[$j] or die "too few values specified for '\$$slVar'\n".
                                       "all variables in linear matrix must have same number of assignment values\n";                    
                fillVariable($slVar, $masterVals[$j]);
            } 
            invokeSlaveMatrix($firstOfSeries, $slaveFile, $slaveVars, $i + $k + 1); # attempt to recurse to next slaveVar to apply combinatorial matrix
        }
    } else { # invoke slave when no more slaveVars are left to assign for this invocation
        invokeSlave($firstOfSeries, $slaveFile); # invocation uses most recent variable assignments
    }
}
sub invokeSlave { # recurse instruction extraction into the slave instructions file
    my ($firstOfSeries, $slaveFile) = @_;
    my $qMasterFile = $currentInstrsFile{$threadLevel}; 
    $threadLevel++; # keep track of slave nesting level    
    my $iFile = findInstrsFile($slaveFile);
    $iFile eq $qMasterFile and die "$qMasterFile cannot invoke itself\n"; 
    my %tmpVariables = %variables;  
    my %tmpENV = %ENV;   
    
    $embeddedTreeFile = $slaveFile;
     
    unless(getFileModuleCache($iFile)){  # execute if file module lookup not already cached
        $threadCounter++; # implied first thread of slave instructions file
        $currentThread{$threadLevel} = $threadCounter; 
        my $currentInstrsDirHold = $currentInstrsDir;  # keep track of the invoking script directory
        $variables{__Q__MASTER__FILE__} = $qMasterFile;                     
        executeInstructions( $iFile ); # act on slave instructions
        setFileModuleCache($iFile);  # maintain a cache of file module lookups
        $currentInstrsDir = $currentInstrsDirHold;  
    }
    my $fos = $$firstOfSeries; # ensure that all return variable use the same first of series value
    foreach my $varRef(@{$returnVariables{$iFile}}){ # retrieve slave return values
        my ($varName, $isReference) = @$varRef;
        defined $variables{$varName} or next; # return variable never actually set by slave
        $$firstOfSeries = 0; # if any variable set, assume all were; otherwise this invocation doesn't change series counter
        if ($fos){ 
            $tmpVariables{$varName} = $variables{$varName}; 
            $tmpENV{$varName} = $ENV{$varName};
        } else { # after the first invocation of a multi-invoke, add this invocation value to a growing list
            if($isReference){  # special behavior when preserving through a referenced variable, i.e. $$varNam
                $tmpVariables{$varName} = $variables{$varName};  # values are NOT appended, will have last encountered value
                $tmpENV{$varName} = $ENV{$varName}; 
            } else {
                $tmpVariables{$varName} = join(" ", $tmpVariables{$varName}, $variables{$varName});   
                $tmpENV{$varName} = join(" ", $tmpENV{$varName}, $ENV{$varName});   
            } 
        }
    }
    delete $returnVariables{$iFile};
    %variables = %tmpVariables; # variables always reset upon return from slave
    %ENV = %tmpENV;
    $currentThread{$threadLevel} = undef; # reset slave level upon return
    $threadLevel--;
}
#------------------------------------------------------------------------
sub setFileModuleCache {  # maintain a cache of file module lookups
    my($iFile) = @_;
    $cacheSlaves{$iFile} or return;
    my $cacheKey = $cacheSlaves{$iFile}[0];
    my $keyValue = $variables{$cacheKey};
    foreach my $varName(@{$cacheSlaves{$iFile}[1]}){
        ${$fileModuleCache{$iFile}{$keyValue}}{$varName} = $variables{$varName};
    }
}
sub getFileModuleCache {  # retrieve from the cache of file module lookups
    my($iFile) = @_;
    $cacheSlaves{$iFile} or return;
    my $cacheKey = $cacheSlaves{$iFile}[0];
    my $keyValue = $variables{$cacheKey};
    my $cache = $fileModuleCache{$iFile}{$variables{$cacheKey}};
    $cache or return;
    foreach my $varName(@{$cacheSlaves{$iFile}[1]}){
        $variables{$varName} = $$cache{$varName}; 
        $ENV{$varName} = $$cache{$varName}; 
        push @{$returnVariables{$iFile}}, [$varName];  
    }   
    return 1;
}
#========================================================================

#========================================================================


# submission-time execution of shell commands
#------------------------------------------------------------------------
sub runImmediately { # when instructed, run simple action commands immediately in q submit shell
    my ($getValue, $shellCommand, @arguments) = @_;
    $shellCommand = join(" ", $shellCommand, @arguments); 
    $shellCommand =~ s/\$NULL//g;    # remove null placeholder
    $shellCommand =~ s/\\(\\*)/$1/g; # remove q's variable escape character    
    if($getValue){
        my $value = qx/$shellCommand/;
        $value =~ s/\n$//;
        $value =~ s/\s+/ /g; 
        return $value;
    } else {
        if($options{'dry-run'}){  # only allow a few essential commands during dry-run
            unless( ($shellCommand =~ m/^echo/ and !$options{'_suppress-echo_'}) or
                     ($shellCommand =~ m/^touch/ or
                      $shellCommand =~ m/^mkdir/) ){
                return;
            }
        }  
        my $result = system($shellCommand);
        $result and die "shell command '$shellCommand' returned an error; no more jobs will be queued\n";
    }
}
#========================================================================

#========================================================================
# process target scripts for queuing
#------------------------------------------------------------------------
sub setQsubMode {  # determine the format of the qsub call, single or matrix
    my ($genericLine, @args) = @_;  # must be done prior to variable replacement
    $qsubMode = 0;
    defined $args[0] or return;
    if($args[0] =~ m|^-|){
        $qsubMode = 1;  # usage: qsub --options scriptFile
    } else {
        $qsubScriptFile = shift @args;
        if($args[0]){
            $qsubMode = 3;  # usage: qsub scriptFile argument [...]
            my @slaveVars;
            foreach my $arg(@args){
                if($arg =~ m|^\$(.+)| and !(defined $variables{$1})){
                    $qsubMode = 4;  # usage: qsub scriptFile invocation_matrix
                    push @slaveVars, " \$$1";
                }
            }
            $qsubMode == 4 and $virtualFiles{$qsubInvokeTarget} = "qsub $qsubScriptFile" . join(" ", @slaveVars);
        } else {
            $qsubMode = 2;  # usage: qsub scriptFile
        }
    }
}
sub parseQsubArguments {  # extract the qsub information based on its format
    my ($genericLine, @args) = @_;
    $qsubMode or die "qsub target script not specified in line:\n$genericLine\n";
    my ($inScript, $qsubOptions, $arguments);
    if($qsubMode == 1){  # usage: qsub --options scriptFile
        $inScript = pop @args;
        $qsubOptions = join(" ", @args); 
        $arguments = "";
    } elsif($qsubMode == 2 or $qsubMode == 3){  # usage: qsub scriptFile [argument ...]
        $inScript = shift @args;
        $qsubOptions = "";
        $arguments = join(" ", @args); 
    } else {
        die "parseQsubArguments error: unknown qsubMode\n";  # qsubMode 4 captured and handled as invoke
    }
    return ($inScript, $qsubOptions, $arguments);
}
#------------------------------------------------------------------------
sub queueScript { # the point of it all - commit target scripts for execution
    my ($currentThread, $genericLine, @args) = @_;
    my ($inScript, $arguments);
    ($inScript, $qsubOptions, $arguments) = parseQsubArguments($genericLine, @args);
    $qsubOptions =~ m/-V/ or  $qsubOptions .= " -V"; # ensure that all environment variables are passed to qsub job
    $qInUse and $qInUse eq 'SGE' and $qsubOptions .= " -terse"; # causes SGE to return only the job ID as output
    %directives = ();
    $jobName = $qsubOptions =~ m/-N\s+(\S+)/ ? $1 : undef; # jobName specified by calling qsub line takes precedence
    $runDir = ($qsubOptions =~ m/-wd\s+(\S+)/ or $qsubOptions =~ m/-d\s+(\S+)/) ? $1 : undef;   
    ($ovrideJobName, $ovrideRunDir, $isJobNameTag, $isRunDirTag) = (defined $jobName, defined $runDir, undef, undef);
    $array = undef;
    my $success = readEnvScripts(\my@outScriptLines, \my@command);
    $success or return;  # environment script was aborted by exitIf or exitUnless
    ($shebang, $isBash) = ($bashPath, 1);
    my $isVirtual = readInScript($inScript, \@outScriptLines, \@command, $genericLine);
    defined $isVirtual or return;  # script was aborted by exitIf or exitUnless    
    $isVirtual or ($currentScriptFile =~ m|.*/(.+)$| and $inScript = $1);  
    $jobName or $jobName = substr($inScript, 0, $qInUse eq 'PBS' ? 15 : 30);  # user never declared a job name
    checkJobName(); 
    $runDir or $runDir = getcwd;    
    my $command = join(" ", @command, $arguments);
    $getOutputFiles and return;  # executing instructions to extract protect or backup fileGlobs; take no action on jobs
    setJobPaths($runDir);     
    checkExtendability($command) or return;    
    $options{'verbose'} and print "queueing job $jobName\n";
    $currentJob++; # jobs are simply numbered in order for tracking    
    updateDependencies($currentThread);  
    finalizeOutScript(\@outScriptLines, $currentScriptFile, $isVirtual, $arguments);
    my $outScript = writeOutScript(\@outScriptLines);    
    addJob($outScript, $command, $arguments);
    $lastJobInThread{$currentThread} = {$currentThread => $currentJob}; # job acts as its own first master for dependency tracking purposes
    updateMasterDependencies($threadLevel, $currentThread); # always update masters when slaves change       
}
sub readEnvScripts { # first read and parse any provided environment scripts
    my ($outScriptLines, $command) = @_;
    foreach my $envFile( getEnvFiles('sh') ){
        -e $envFile and (defined readInScript($envFile, $outScriptLines, $command, undef, 1) or return);   
    }   
    return 1;  
}
sub readInScript { # then read and parse target script
    my ($inScript, $outScriptLines, $command, $genericLine, $isEnv) = @_;
    $inScript = getInscript($inScript, $genericLine);   
    my $inScriptLines;
    if(ref($inScript) eq 'SCALAR'){
        $inScriptLines = slurpFile($inScript); 
    } elsif($fileContents{$inScript}){
        $inScriptLines = $fileContents{$inScript};
    } elsif(-e $inScript){
        $inScriptLines = slurpFile($inScript); 
        $fileContents{$inScript} = $inScriptLines;  # maintain a cache of encountered script files
    } else {
        return;
    }
    open my $inH, "<", \$inScriptLines or die "could not open $inScript for reading: $!\n";
    while (my $scriptLine = <$inH>){
        chomp $scriptLine;
        $scriptLine =~ s/\r//g;
        if($isEnv){ # environment scripts only pass non-blank non-comment lines (directives are passed)
            $scriptLine =~ m/\S/ or next;
            if($scriptLine =~ m/^#/){ $scriptLine =~ m/^#(\!|q|\$|PBS)/ or next};
        }  
        replaceDirectiveVars(\$scriptLine, $inScript, $scriptLine);
        processScriptLine(\$scriptLine, $inScript, $command) or return;
        if(defined $scriptLine){ # collect script line for output script and tracked composite command
            push @$outScriptLines, $scriptLine;        
            compressScriptLine($scriptLine, $command); 
        }
    }     
    close $inH; 
    return ref($inScript);
}
sub getInscript {
    my ($inScript, $genericLine) = @_;
    $inScript or die "target script not specified in qsub line:\n$genericLine\n"; 
    if($virtualFiles{$inScript}){
        $currentScriptFile = '';
        return \$virtualFiles{$inScript}; # virtual files declared in master (or slave) take precedence over disk files 
    } 
    my @dirs =  ($masterDir,         # look first in same directory as the master instructions file
                 $currentInstrsDir,  # then look in the same directory as the invoking instructions file
                 $modulesDir,        # then look in the q modules directory
                 '',                 # then attempt to interpret as relative or absolute file path
                 );  
    foreach my $i(0..$#dirs){ 
        defined $dirs[$i] or next; 
        my $file = $dirs[$i] ? "$dirs[$i]/$inScript" : $inScript;  
        if($fileContents{$file} or -e $file){
            $file =~ m|^\/| or $file = abs_path($file); 
            $currentScriptFile = $file;
            return $file;
        }
    }    
    die "could not find script $inScript specified in qsub line:\n$genericLine\n";     
}
sub replaceDirectiveVars { # allow job scheduler directives to use variable values inherited from q
    my ($scriptLine, $inScript, $genericLine) = @_;
    if($$scriptLine =~ m/^#(\$|PBS)\s+\S+/){
        while($$scriptLine =~ m/\$(\w+)/){ 
            my $varName = $1; 
            my $varValue = $variables{$varName};
            defined $varValue or die "variable \$$varName used before being assigned a value in $inScript line:\n$genericLine\n";
            $varValue =~ s/\$NULL//g;    # remove null placeholder
            $$scriptLine =~ s/\$$varName/$varValue/g;
        }
        $$scriptLine =~ s/\\(\\*)/$1/g; # remove q's variable escape character
    }
}
sub processScriptLine { # act on directives in environment and target scripts
    my ($scriptLine, $inScript, $command) = @_;
    if($$scriptLine =~ m/^#q\s+require\s+(.+)/){ # line specifies q required variable(s)
        my $requiredVars = $1;
        $requiredVars = checkRequiredVars($inScript, $requiredVars, $command);
        push @{$directives{q}}, "require $requiredVars";
        $$scriptLine = undef;
    } elsif($$scriptLine =~ m/^#q\s+option\s+(.+)/){ # line specifies q optional variable(s)
        my $optionalVars = $1;
        $optionalVars = checkOptionalVars($inScript, $optionalVars, $command);
        push @{$directives{q}}, "option $optionalVars";
        $$scriptLine = undef;      
    } elsif($$scriptLine =~ m/^#q\s+(exitIf\s+.+)/ or $$scriptLine =~ m/^#q\s+(exitUnless\s+.+)/){ # allow abort out of script file
        processInstruction($1, $inScript) or return;      
    } elsif($$scriptLine =~ m/^#q\s+.+/){ # fail on all other attempted q directives
        die "unknown q directive in line:\n$$scriptLine\n";
    } elsif($$scriptLine =~ m/^(#\$\s+-wd\s+)(.+)/ or $$scriptLine =~ m/^(#PBS\s+-d\s+)(.+)/){ # collect requested execution Directory
        if($ovrideRunDir){ $$scriptLine = "$1$runDir" } else { $runDir = $2 }
        $isRunDirTag = 1; 
    } elsif($$scriptLine =~ m/^(#\$\s+-N\s+)(.+)/ or $$scriptLine =~ m/^(#PBS\s+-N\s+)(.+)/){ # collect job name specified in script
        if($ovrideJobName){ $$scriptLine = "$1$jobName" } else { $jobName = $2 }
        $isJobNameTag = 1;
    } elsif($$scriptLine =~ m/^#\$\s+-t\s+(.+)/ or $$scriptLine =~ m/^#PBS\s+-t\s+(.+)/){ # determine if this is an array job
        #$array = 1;  # versions prior to 0.4.0 used simple boolean to indicate an array job
        $array = getArrayTaskIDs($1, $scriptLine);  # later version collect the complete set of array task IDs
    } elsif($$scriptLine =~ m/^#!(.+)/){ # manage qsub target shebang line
        $shebang = $1;  
        -x $shebang or die "unrecognized program target in shebang line:\n$$scriptLine";
        if($shebang =~ m/bash$/){ $directives{'!'} = $shebang } else { $isBash = undef }
        $$scriptLine = undef;
    } 
    my %disallowedOptions = map { $_ => 1 } qw(j o e);  # output path directives must be set by q! 
    if(defined $$scriptLine and $$scriptLine =~ m/^#(\$|PBS)\s+-(\S+)/){ 
        my ($qTag, $option) = ($1, $2);  # keep track of directives to allow them to be nicely ordered at the top of the script
        if($disallowedOptions{$option}){
            push @{$directives{$qTag}}, "# directive overridden by q:  $$scriptLine";
        } else {
            push @{$directives{$qTag}}, $$scriptLine;
        }   
        $$scriptLine = undef;  
    }     
    return 1;                                   
}
sub getArrayTaskIDs {  # determine the task IDs that will be generated by an array job
    my ($tOption, $scriptLine) = @_;
    $getOutputFiles and return 1;
    my $usage = "invalid array specification: $$scriptLine\n";
    my ($start, $end, $step) = (1, 1, 1);
    if($tOption =~ m/^(\d+)$/){
        $start = $1;
        $end = $start;
    } elsif($tOption =~ m/^(\d+)-(\d+)$/){
        ($start, $end) = ($1, $2);
    } elsif($tOption =~ m/^(\d+)-(\d+):(\d+)$/){
        ($start, $end, $step) = ($1, $2, $3);
    } else {
        die $usage;
    }
    $start >= 1 or die $usage;
    $end >= 1 or die $usage;
    $end >= $start or die $usage;
    my @taskIDs;
    for (my $i = $start; $i <= $end; $i += $step){ push @taskIDs, $i }
    return join(",", @taskIDs);  # task IDs stored as comma-delimited list
}    

sub checkRequiredVars { # when user provides variable syntax requirements, make sure called q files have specified them
    my ($inScript, $requiredVars, $command) = @_; 
    $requiredVars =~ s/,//g;
    $requiredVars =~ s/\s//g;
    my @requiredVars = split('\$', $requiredVars);
    shift @requiredVars; 
    $requiredVars = "\n";
    foreach my $varName(@requiredVars){
        my $varValue = $variables{$varName};
        defined $varValue or die ref($inScript)?"embedded script":$inScript," requires a value for variable \$$varName\n"; 
        $varValue =~ s/\$NULL//g;    # remove null placeholder
        $varValue =~ s/\\(\\*)/$1/g; # remove q's variable escape character"";
        my $scriptLine = "$varName=\"$varValue\"";
        $requiredVars .= "$scriptLine\n";   
        push @$command, $scriptLine;   
    }
    return $requiredVars;
}
sub checkOptionalVars { # set user specific values for optional script variables, otherwise set to default
    my ($inScript, $optionalVars, $command) = @_; 
    my $outputVars = "\n";    
    while($optionalVars =~ s|\s*\$(\w+)\s*\[(.*?)\]||){
        my ($varName, $defaultValue) = ($1, $2);
        defined $defaultValue or $defaultValue = "";
        my $varValue = $variables{$varName};
        defined $varValue or $varValue = $defaultValue;
        $varValue =~ s/\$NULL//g;    # remove null placeholder
        $varValue =~ s/\\(\\*)/$1/g; # remove q's variable escape character"";        
        my $scriptLine = "$varName=\"$varValue\"";
        $outputVars .= "$scriptLine\n";   
        push @$command, $scriptLine;           
    }
    $optionalVars =~ s|\s||g;
    $optionalVars and die ref($inScript)?"embedded script":$inScript," q option directive is malformed\n";  
    return $outputVars;
}
sub compressScriptLine { # minimize script line to its essential elements for command tracking
    my ($scriptLine, $command) = @_;
    $scriptLine =~ s/#.*$//; # compress line to minimal command string - ignore comments    
    $scriptLine =~ m/\S/ or return;
    $scriptLine =~ s/\n/ /g; # collapse to single-space delimited word series
    $scriptLine =~ s/\t/ /g; # collapse to single-space delimited word series
    $scriptLine =~ tr/ / /s;      
    $scriptLine =~ s/\s+$//; # strip trailing white space
    $scriptLine =~ s/^\s+//; # strip leading white space
    push @$command, $scriptLine;      
}
sub finalizeOutScript { # assemble the final script for this job using information from qsub options and environment and target scripts
    my ($outScriptLines, $inScriptPath, $isVirtual, $arguments) = @_; # assemble in reverse so that commands are stacked on top of use script lines
    processEmbeddedTarget($outScriptLines, $inScriptPath, $isVirtual, $arguments); 
    addDirectives('PBS', 'd', 'N', 'W', undef, 'depend=afterok:'); # put together and organize the complete set of scheduler directives
    addDirectives('$', 'wd', 'N', 'hold_jid', 1, ''); 
    $directives = assembleDirectives('PBS')."\n";
    $directives .= assembleDirectives('$');
    unshift @$outScriptLines, $directives; # place one copy of the directives in target script for future reference (these aren't actually read by qsub)
    $directives{'q'} and unshift @$outScriptLines, "#q ".join("\n#q ", @{$directives{q}}); # q directives
    unshift @$outScriptLines, "getTaskID\n";                         # collect TASK_ID on user's behalf by default; must occur before q requires    
    unshift @$outScriptLines, "checkPredecessors";                   # check for timed out predecessors on users behalf 
    unshift @$outScriptLines, "source \"$qDir/lib/utilities.sh\""; # make q script utilities available to target scripts   
    unshift @$outScriptLines, "echo";
    unshift @$outScriptLines, "echo \"q: running on host: \$HOSTNAME\"";
    $directives{'!'} or $directives{'!'} = $bashPath;  
    unshift @$outScriptLines, "#!".$directives{'!'}."\n"; # finish with shebang line  
}
sub processEmbeddedTarget{  # replace script contents with a call to non-bash scripts
    my ($outScriptLines, $inScriptPath, $isVirtual, $arguments) = @_;
    $isBash and return;  # embedded bash targets are simply incorporated into q-parsed script
    if($isVirtual){  # create a hidden disk copy of embedded non-bash script targets
        abs_path($currentInstrsFile{$threadLevel}) =~ m|(.+)/(.+)| or die "finalizeOutScript: error parsing current instrsFile name\n";
        my ($qPath, $qFile) = ($1, $2);
        $inScriptPath = "$qPath/.q.embedded";
        mkdir $inScriptPath;
        $inScriptPath = "$inScriptPath/$qFile.$qsubScriptFile";
        open my $outH, ">", $inScriptPath or die "could not open $inScriptPath for writing: $!\n";
        print $outH join("\n", @$outScriptLines);
        close $outH;    
    }
    @$outScriptLines = "$shebang $inScriptPath $arguments\n";
}
sub addDirectives { # include directives established by q = runDir, jobName, dependencies
    my ($qTag, $runDirOption, $jobNameOption, $dependOption, $dependSubsDiv, $dependPrefix) = @_; 
    $isRunDirTag or  push @{$directives{$qTag}}, "#$qTag -$runDirOption $runDir";  
    $isJobNameTag or push @{$directives{$qTag}}, "#$qTag -$jobNameOption $jobName";  
    if($dependencies){    
        $dependSubsDiv and $dependencies =~ s/:/,/g;
        push @{$directives{$qTag}}, "#$qTag -$dependOption $dependPrefix$dependencies";  
    }  
}
sub assembleDirectives { # create organized blocks of directives for SGE and PBS
    my ($qTag) = @_;
    $directives{$qTag} or return "";
    return join("\n", @{$directives{$qTag}})."\n";
}
sub writeOutScript { # put it all together and commit the script for this job
    my ($outScriptLines) = @_;
    $options{'dry-run'} and return;
    my $scriptDir = getScriptDir($qInUse);
    my $scriptBase = "$scriptDir/$jobName";      
    $scriptBase .= ".sh";   
    my $outScript = $scriptBase;    
    my $i = 0;    
    while (-e $outScript){
        $i++;
        $outScript = "$scriptBase.$i";
    }
    my $finalOutScript = join("\n", @$outScriptLines);
    $finalOutScript =~ s/\n\s+\n/\n\n/g; # remove extra blank lines for pretty parsing
    open my $outH, ">", $outScript or die "could not open $outScript for writing: $!\n";  
    print $outH $finalOutScript;
    close $outH;   
    return $outScript;
}
sub setJobPaths { # override any output paths specific by user; q needs to know where things are!
    my ($runDir) = @_;
    $runDir or $runDir = getcwd;   
    -d $runDir or qx/mkdir -p $runDir/;
    my $logDir = getLogDir($qInUse);
    $ENV{Q_LOG_DIR} = $logDir;
    push @{$directives{'$'}}, "#\$    -j y"; 
    push @{$directives{'PBS'}}, "#PBS  -j oe"; 
    push @{$directives{'$'}}, "#\$    -o $logDir"; 
    push @{$directives{'PBS'}}, "#PBS  -o $logDir";      
} 
sub checkJobName {
    $jobName =~ m|\s| and die "job name cannot contain white space: '$jobName'\n";
    #$qInUse and $qInUse eq 'PBS' and length($jobName) > 15 and die "PBS job names cannot exceed 15 characters:  '$jobName'\n";
} 
sub addJob { # act on the assembled job
    my ($targetScript, $command, $arguments) = @_;
    my $jobID = submitJob($targetScript, $arguments);
    $options{'dry-run'} or saveEnvironment($qInUse, $jobID);
    $jobIDs{$currentJob} = $jobID;
    push @statusInfo, [$jobName, $currentJob, $targetScript, $command, $array, $currentInstrsFile{$threadLevel}, $currentScriptFile];
    unless($options{'_suppress-echo_'} or $qInUse eq 'local'){
        my $pred = join(",", @{$jobPred{$currentJob}});
        $pred or $pred = "";
        my $jobID = $options{'dry-run'} ? 0 : $jobIDs{$currentJob}; 
        padSubmitEchoColumns($jobName, $array ? '@' : ' ', $jobID, $currentJob, $pred);   
    }
    return $jobID;
}
sub padSubmitEchoColumns {  # ensure pretty parsing of echoed status table
    my(@in) = @_;
    my @columnWidths = (30, 1, 7, 5, 20);
    foreach my $i(0..4){
        my $value = $in[$i];
        my $outWidth = $columnWidths[$i];    
        $value =~ s/\s+$//;
        $value = substr($value, 0, $outWidth);        
        my $padChar = " ";
        $i == 0 and $value .= " " and $padChar = "-"; 
        my $inWidth = length($value);
        $inWidth < $outWidth and $value .= ($padChar x ($outWidth - $inWidth ));
        print "$value"."  ";
    }
    print "\n";
}
#========================================================================

#========================================================================
# job and thread dependencies
#------------------------------------------------------------------------
sub updateDependencies { # assemble the current dependency chain for the job about to be queued
    my ($thread) = @_;
    $qInUse eq 'local' and return; # local is just a synchronous series of jobs, dependencies are irrelevant
    my @pred = getLastJobsInThread($thread); # jobs already exist in this thread to depend on
    scalar(@pred) or ($threadPred{$thread} and @pred = @{$threadPred{$thread}}); #otherwise obey any dependencies of the current thread
    my $i = $threadLevel - 1; # master thread to the current slave
    until(scalar(@pred) or $i < 0){ # otherwise depend on last jobs of an invoking thread
        my $masterThread = getCurrentThread($i);
        @pred = getLastJobsInThread($masterThread, 1); # special call to held invoke dependencies; last job in master thread
        scalar(@pred) or ($threadPred{$masterThread} and @pred = @{$threadPred{$masterThread}}); # otherwise obey dependencies of the master thread
        $i--; # recurse upwards through thread levels in case invoking slave thread has declared no prior jobs
    }
    @{$jobPred{$currentJob}} = sort {$a <=> $b} @pred; # order and hold jobs for screen echo
    my %pred; # use hash to assemble final dependency list to prevent redundant jobs in list
    foreach my $job(@pred){
        $pred{$jobIDs{$job}}++; #convert predecessor job numbers to job IDs
        $jobSucc{$job}{$currentJob}++; #track successors as reciprocal of predecessors
    }
    if($options{depend}){ # add any jobIDs supplied by the user
        foreach my $dependID(split(",", $options{depend})){ $pred{$dependID}++ }
    }
    $dependencies = join(':',  sort {$a <=> $b} keys %pred); # return a numerically sorted, colon-delimited list of unique jobIDs
    $ENV{Q_PREDECESSORS} = $dependencies;
}
sub updateMasterDependencies { # whenever a slave dependency is modified, recurse the change up the dependency chain to all masters/ancestors
    my ($slaveLevel, $slaveThread) = @_;
    if ($slaveLevel > 0){ # is a slave thread
        updateMasterDependencies_($slaveLevel - 1, $slaveThread);
        my $masterThread = getCurrentThread($slaveLevel - 1);

        foreach my $mergeThread(@{$mergeThreads{$slaveThread}}){   # master no longer need to depend on merging threads
            delete $lastJobInThread{$masterThread}{$mergeThread};  # once jobs are declared in a continuation thread
            updateMasterDependencies_($threadLevel - 2, $masterThread);
        }
    }
}
sub updateMasterDependencies_ { # replace dependency information in master of modified slave
    my ($masterThreadLevel, $slaveThread) = @_;
    while($masterThreadLevel >= 0){
        my $masterThread = getCurrentThread($masterThreadLevel);
        %{$lastJobInThread{$masterThread}{$slaveThread}} = %{$lastJobInThread{$slaveThread}}; 
        $slaveThread = $masterThread;
        $masterThreadLevel--; # recurse up the master/ancestor list
    }
}
#========================================================================

#========================================================================
# submit job to queue
#------------------------------------------------------------------------
sub submitJob{ # disperse the job as indicated by qInUse
    my ($targetScript, $arguments) = @_;
    $options{'dry-run'} or qx/chmod u+x $targetScript/;
    if ($qInUse eq 'local'){   
        my @taskIDs = $array ? split(",", $array) : ("1");
        my $jobID;
        foreach my $taskID(@taskIDs){  # array jobs are submitted serially when executed locally
            $ENV{PBS_ARRAYID} = $taskID;  # ensure that taskID is passed to job
            $ENV{SGE_TASK_ID} = $taskID;
            $jobID = submitLocal($targetScript, $arguments);  # just return the last local job id in the array
            $ENV{PBS_ARRAYID} = undef;
            $ENV{SGE_TASK_ID} = undef;
        }
        return $jobID;
        #$array and die "array jobs are not compatible with option --execute\n";
        #return submitLocal($targetScript, $arguments);  
    } else {
        return submitQueue($targetScript, $arguments);
    }
}
sub getExecutionFile {      # wrap user script in q-defined execution script that uses GNU time utility
    my ($targetScript, $arguments) = @_;  # to collect, parse and format exit_status and resource utilization   
    my $executionFile = "$masterFile.executing.sh"; # create temporary script file
    open my $outH, ">", $executionFile or die "could not open $executionFile for writing: $!\n";
    my $date = '`date +\'%a %D %R\'`';
    print $outH '#!', $bashPath, "\n";
    print $outH $directives; #add the second instance of directives to execution script - these are the ones processed by qsub
    print $outH 'echo "q: target script: ', $targetScript, '"', "\n";
    print $outH 'echo "q: execution started: ', $date, '"', "\n";
    print $outH $timePath, ' -f "\nq: exit_status: %x; walltime: %E; seconds: %e; maxvmem: %MK; swaps: %W; I/O: %I/%O;" ', 
                $targetScript, " ", $arguments, "\n";
    print $outH 'x=$?', "\n";
    print $outH 'echo "', $memoryMessage, '"', "\n";
    print $outH 'echo "q: execution ended: ', $date, '"', "\n";
    print $outH '[ "$x" -gt 0 ] && x=100', "\n";    
    print $outH 'exit $x', "\n"; # force all q exit statuses to be 0 or 100 for appropriate SGE management
    close $outH; 
    return $executionFile; 
}
sub submitLocal { # run the script in shell if queue is suppressed
    my ($targetScript, $arguments) = @_;
    my $separatorLength = $options{'dry-run'} ? 0 : 75;
    $options{'dry-run'} or print "=" x $separatorLength, "\n";
    $options{'_suppress-echo_'} or print "$jobName\n"; 
    $options{'dry-run'} or print "~" x $separatorLength, "\n";
    my $jobID = 0;
    unless($options{'dry-run'}){
        my $executionFile = getExecutionFile($targetScript, $arguments);
        qx/chmod u+x $executionFile/;
        my $logContents = qx/$executionFile 2>&1/; # execute temporary script file, merging stderr to stdout
        unlink $executionFile;         
        print $logContents;   
        $jobID = getLocalJobID();
        my $jName = $jobName;
        $jName =~ s/\s+$//;
        my $logFile = getLogDir('local')."/$jName.o$jobID"; 
        open my $logFileH, ">", $logFile or die "could not open $logFile for writing: $!\n"; 
        print $logFileH "$logContents\n";
        close $logFileH;
        $jobInfos{$currentJob} = {};
        parseLogFile($jobInfos{$currentJob}, $logContents, 1);
        (defined $jobInfos{$currentJob}{exit_status} and $jobInfos{$currentJob}{exit_status} == 0) 
            or die "=" x $separatorLength."\n\njob error: no more jobs will be queued\n";    
        $jobsAdded = 1;          
    }
    $options{'dry-run'} or print "=" x $separatorLength, "\n";
    return $jobID; 
}
sub getLocalJobID {
    my $localDir = getLogDir('local');
    my @logFiles = <$localDir/*.o*>;
    my $maxJobID = 0;
    foreach my $logFile(@logFiles){
        $logFile =~ m|$localDir/.+\.o(\d+)| or next;
        $maxJobID >= $1 or $maxJobID = $1;
    }
    return $maxJobID + 1;
}
sub submitQueue { # submit to qsub
    my ($targetScript, $arguments) = @_;
    $options{'dry-run'} and return $currentJob; # dry-run just returns q internal job number
    my $executionFile = getExecutionFile($targetScript, $arguments);  
    my $qsubCommand = "$schedulerDir/qsub $qsubOptions $executionFile";
    my $jobID = qx/$qsubCommand/;
    $jobID =~ m/^(\d+).*/ or die "error recovering qsub jobID\n";
    $jobID = $1;
    $jobsAdded = 1;
    unlink $executionFile;
    return $jobID;  
}
#========================================================================

#========================================================================
# generate status file
#------------------------------------------------------------------------
sub generateStatusFile { # create the file that is used by subsequent q commands such as delete, etc.
    my $time = getTime();
    my $createdArchive = archiveStatusFiles();  # attempt to create an archive copy of a pre-existing status file 
    open my $statusFileH, ">>", $statusFile or die "could not open $statusFile for appending: $!\n";
    my $user = $ENV{USER};
    print $statusFileH 
        "qType\t$qInUse\n",
        "submitted\t$time\t$user\n",
        join("\t", qw(  jobName
                        jobID
                        array
                        jobNo
                        predecessors
                        successors
                        start_time
                        exit_status
                        walltime
                        maxvmem   
                        targetScript
                        command
                        instrsFile
                        scriptFile
                        user
                        qType ))."\n";
    foreach my $jobInfo(@statusInfo){
        my ($jobName, $job, $targetScript, $command, $array, $instrsFile, $scriptfile) = @$jobInfo;
        my ($startTime, $exitStatus, $wallTime, $maxVmem) = ('', '', '', '');
        if($jobInfos{$job}){ # job executed locally; already have status information
            $startTime = $jobInfos{$job}{start_time};
            $exitStatus = $jobInfos{$job}{exit_status};
            $wallTime = $jobInfos{$job}{walltime};
            $maxVmem = $jobInfos{$job}{maxvmem};   
        }
        $array or $array = '';
        my $pred = $jobPred{$job} ? join(",", @{$jobPred{$job}}) : '';
        my $succ = $jobSucc{$job} ? join(",", sort {$a <=> $b} keys %{$jobSucc{$job}}) : '';
        my $jobID = $jobIDs{$job};
        print $statusFileH 
            join("\t",  $jobName,
                        $jobID,
                        $array,
                        $job,
                        $pred,
                        $succ,
                        $startTime, 
                        $exitStatus, 
                        $wallTime, 
                        $maxVmem,
                        $targetScript,
                        $command,
                        $instrsFile,
                        $scriptfile,
                        $user,
                        $qInUse )."\n";
    }
    close $statusFileH;
    $createdArchive or archiveStatusFiles();  # for 1st write of status file, create an immediate archive of it
}
#========================================================================

#========================================================================
# pipeline management automatic jobs
#------------------------------------------------------------------------
sub setAutoUnprotect {  # communicate whether the pipeline will protect files that need unprotecting to allow proper execution
    $options{'_suppress-echo_'} or return;  # only set unprotect marker at the end of a syntax check for an execution run
    my $autoUnprotectFile = getAutoUnprotectFile();     # marker is passed as the existence of a file
    $needAutoProtect and isProtectedOutput() and qx|touch $autoUnprotectFile|;  # no other mechanism for passing from one q instance to another
}
sub needAutoUnprotect {  # determine whether q unprotect is needed prior to queuing other jobs
    $options{'dry-run'} and return 0;  # no action needed for dry run
    my $autoUnprotectFile = getAutoUnprotectFile();
    -e $autoUnprotectFile or return 0;  # no protect commands were encountered during the syntax check
    unlink $autoUnprotectFile; 
    getPermission("$masterFileName protects file(s) that already exist and are write-protected\n".
                  "continuing will prepend a 'q unprotect' job to allow output files to be re-written", 1) and return 1;
    die "aborting\n";
}
sub getAutoUnprotectFile {
    return "$masterDir/.$masterFileName.autoUnprotect";
}
sub addAutoFileGlob {  # collect the files specified by a protect command, etc.
    my ($genericLine, $filesHash, @line) = @_;
    $line[0] or die "no file glob specified in line: $genericLine\n";
    if($line[0] eq 'find'){  # allow use of 'find' statement to specify file list
        my $find = join(" ", @line);
        $find =~ s/\\(\\*)/$1/g; # remove q's variable escape character   
        $$filesHash{$find}++;
    } else {  # otherwise assume a simple fileGlob, or list of fileGlobs, was provided
        foreach my $fileGlob(@line){
            $fileGlob =~ s/\\(\\*)/$1/g; # remove q's variable escape character   
            $$filesHash{$fileGlob}++;
        }
    }      
}    
sub addAutoDirectory {  # collect the directories specified by a backup command, etc.
    my ($genericLine, $filesHash, @line) = @_;
    $line[0] or die "no directory specified in line: $genericLine\n";
    foreach my $directory(@line){
        $directory =~ s/\\(\\*)/$1/g; # remove q's variable escape character   
        $$filesHash{$directory}++;
    }
}
sub initializeAutoThread {  # handle the hidden parsing of automatic threads present in all job sequences
    my ($timing) = @_;
    my $threadName = "q_auto_$timing";
    my @mergeThreads = keys %masterThreadNumbers;
    initializeThread($masterFile, \%masterThreadNumbers, 'thread', $threadName, @mergeThreads);
}     
sub getPassedAutoOptions {  # pass relevant options on to requested automatic jobs  
    my ($command) = @_;
    my $autoOptions = "";
    foreach my $option(keys %{$commandOptions{$command}}){
        defined $options{$option} and $autoOptions .= " --$option $options{$option}";
    }   
    return $autoOptions;
}    
sub addAutoJob {  # add requested automatic jobs before and after user-requested jobs
    my ($needed, $autoType, $mbRAM, $wallTime, $qOptions) = @_;
    $needed or return;   
    my $jobName = "$autoType\_$masterFileName";
    my $scriptFile = "$masterDir/$jobName.sh";         
    $qOptions or $qOptions = "";    
    open my $scriptH, ">", $scriptFile or die "error opening $scriptFile for writing: $!\n";
    print $scriptH join("\n",
        '#!'.$bashPath,
        '#$    -N  '.$jobName,
        '#$    -wd '.$masterDir,
        '#$    -l  vf='.$mbRAM.'M',
        '#$    -l  h_rt='.$wallTime,
        '#PBS  -N  '.$jobName,
        '#PBS  -d  '.$masterDir,
        '#PBS  -l  mem='.$mbRAM.'mb',
        '#PBS  -l  walltime='.$wallTime,        
        "perl $qTarget $autoType $qOptions $masterFile\n");        
    close $scriptH;
    $qsubMode = 2;
    $currentInstrsFile{0} or $currentInstrsFile{0} = $masterFile;
    queueScript(getCurrentThread(0), "qsub $scriptFile", $scriptFile);
    $qsubMode = 0;
    unlink $scriptFile;
}
#========================================================================

1;


