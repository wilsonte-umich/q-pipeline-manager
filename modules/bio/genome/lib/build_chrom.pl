#!/usr/bin/perl
use strict;
use warnings;
use Cwd(qw(abs_path));

#q    require $INPUT_GENOME_DIR $INPUT_GENOME $LIFT_BEDS $CHANGE_FORMAT
#q    require $OUTPUT_GENOME_DIR $OUTPUT_GENOME $CHANGES_FILE $CHROMOSOME $FEATURES_FILE

#$    -N  build.$OUTPUT_GENOME.$CHROMOSOME
#$    -wd $GENOMES_DIR
#$    -l  vf=4G
#$    -l  h_rt=5:00:00

#PBS  -N  build.$OUTPUT_GENOME.$CHROMOSOME
#PBS  -d  $GENOMES_DIR
#PBS  -l  mem=4gb
#PBS  -l  walltime=5:00:00

# define variables
my $inGenome = $ENV{INPUT_GENOME};
my $outGenome = $ENV{OUTPUT_GENOME};
my $chrom = $ENV{CHROMOSOME};
print "building $outGenome version of $chrom from genome $inGenome\n";
my $inDir = $ENV{INPUT_GENOME_DIR};
my $outDir = $ENV{OUTPUT_GENOME_DIR};
my $inFasta =  "$inDir/chromosomes/$chrom.fa";
my $outFasta = "$outDir/chromosomes/$chrom.fa";
-e $inFasta or die "build_chrom.pl error: file not found: $inFasta\n";
my $changeFormat = $ENV{CHANGE_FORMAT};
($changeFormat eq 'COORDINATES' or $changeFormat eq 'FLANKS') or die "build_chrom.pl error: unkown change format: $changeFormat\n";
my $changeSub = $changeFormat eq 'COORDINATES' ? \&COORDINATES : \&FLANKS;
my $changes = "slurp $ENV{CHANGES_FILE} | awk '\$1==\"$chrom\"'";
my $sortedChanges = $changeFormat eq 'COORDINATES' ? "$changes | sort -k2,2nr" : $changes;
my $chromSeq;
my %inLiftBeds = map { $_ => 1 } split(" ", $ENV{LIFT_BEDS});
my @inLiftBeds = keys %inLiftBeds;
my ($liftStartI, $liftEndI, %liftStartIs, %liftEndIs, %revLiftStartIs, %revLiftEndIs, 
    @liftStartIs, @liftEndIs, %liftFeatures, %liftCombIs, %outBedHs) = (-1, -1);
my $features = "slurp $ENV{FEATURES_FILE} | awk '\$2==\"$chrom\"'";

# main execution block
if(isStreamData($changes)){
    processChanges();
    printRevisedChrom();
    printLiftBeds();
} else {
    printAsIs();
}
if(isStreamData($features)){
    addNewFeatures();
} else {
    print "no new $chrom features requested\n";
}
print "done\n";

# check for stream data
sub isStreamData {
    my ($stream) = @_;
    my $isStreamData = qx/$stream | head -n1 | wc -l/;
    chomp $isStreamData;
    return $isStreamData;
}

# process and print sequence changes
sub processChanges {
    print "changes requested in $chrom, processing\n";
    loadChrom();
    print length($chromSeq)."\t$chrom beginning length\n";
    loadLiftBeds();
    open my $inH, "-|", $sortedChanges or die "build_chrom.pl error: could not open change stream for reading: $!\n";
    while(my $change = <$inH>){   
        chomp $change;  
        $change =~ s/\r//g;
        $change or next;
        $change =~ m/^\s*#/ and next;
        my ($chrom, $before, $after, $toInsert) = split(/\s+/, $change);
        defined $after or next;    
        $toInsert or $toInsert = "";
        $toInsert = uc($toInsert);        
        &$changeSub(\$before, \$after, \$toInsert, length($toInsert), \$change);
    }
    close $inH;
    print length($chromSeq)."\t$chrom ending length\n";
}
sub printRevisedChrom {
    print "printing revised $chrom sequence\n";
    open my $outH, ">", $outFasta or die "build_chrom.pl error: could not open $outFasta for writing: $!\n";
    print $outH ">$chrom\n";
    print $outH join("\n", unpack("(A50)*", $chromSeq)), "\n";
    close $outH;
}
sub printLiftBeds {
    print "printing revised lift BED files for $chrom\n";
    if($liftStartI >= 0){
        foreach my $ij (keys %liftFeatures) {
            my ($startI, $endI) = split(" ", $ij);
            my $start = $revLiftStartIs{$startI};
            $start > 0 or next;            
            my $end = $revLiftEndIs{$endI};
            $end > 0 or next;            
            foreach my $feature(@{$liftFeatures{$ij}}){
                my ($outLiftBed, @trailing) = @$feature;
                my $outH = getOutBedH($outLiftBed);
                print $outH join("\t", $chrom, $start - 1, $end, @trailing), "\n";
            }
        }
    }
    closeOutBedHs();
}
sub printAsIs {
    print "no changes in $chrom, copying files as is\n";
    qx/cp $inFasta $outFasta/;
    foreach my $inLiftBed_(@inLiftBeds){
        my ($inLiftBed, $outLiftBed) = getOutLiftBed($inLiftBed_);
        qx/slurp $inLiftBed | awk '\$1=="$chrom"' > $outLiftBed/;
    }  
}

# load a chromosome for modification
sub loadChrom {
    open my $inH, "<", $inFasta or die "build_chrom.pl error: could not open $inFasta for reading: $!\n";
    while(my $line = <$inH>){
        chomp $line;
        $line =~ s/\s//g;
        $line or next;
        $line =~ m/^>/ and next;
        $chromSeq .= uc($line);
    }
    close $inH;
}
sub loadLiftBeds {
    foreach my $inLiftBed_(@inLiftBeds){
        my ($inLiftBed, $outLiftBed) = getOutLiftBed($inLiftBed_);
        my $features = "slurp $inLiftBed | awk '\$1==\"$chrom\"'";
        open my $inH, "-|", $features or die "build_chrom.pl error: could not open $inLiftBed stream for reading: $!\n";
        while(my $line = <$inH>){
            chomp $line;
            $line =~ s/\r//g;
            $line or next;
            my ($chrom, $start, $end, @trailing) = split("\t", $line);
            $start++;  # this file works in 1-referenced coordinates
            unless(defined $liftStartIs{$start}){
                $liftStartI++;            
                $liftStartIs{$start} = $liftStartI;
            }
            unless(defined $liftEndIs{$end}){
                $liftEndI++;            
                $liftEndIs{$end} = $liftEndI;
            }
            my $ij = "$liftStartIs{$start} $liftEndIs{$end}";
            push @{$liftFeatures{$ij}}, [$outLiftBed, @trailing];
        }
        close $inH;  
    }
    if($liftStartI >= 0){
        %revLiftStartIs = reverse %liftStartIs; 
        %revLiftEndIs = reverse %liftEndIs; 
        @liftStartIs = sort {$revLiftStartIs{$b} <=> $revLiftStartIs{$a}} 0..$liftStartI;
        @liftEndIs = sort {$revLiftEndIs{$b} <=> $revLiftEndIs{$a}} 0..$liftEndI;    
    }
}
sub getOutLiftBed {
    my ($inLiftBed) = @_;
    -e $inLiftBed or $inLiftBed = "$inDir/$inLiftBed";
    my $outLiftBed = $inLiftBed;
    $outLiftBed =~ s/$inGenome/$outGenome/g;
    $outLiftBed .= ".chrom_$chrom";
    return ($inLiftBed, $outLiftBed);
}

# makes the sequence changes
sub COORDINATES {
    my ($lastBefore, $firstAfter, $toInsert, $insertionLength) = @_;
    $chromSeq = join("", substr($chromSeq, 0, $$lastBefore), 
                         $$toInsert, 
                         substr($chromSeq, $$firstAfter - 1));
    liftOver($$lastBefore, $$firstAfter, $insertionLength);    
}
sub FLANKS {
    my ($leftFlank, $rightFlank, $toInsert, $insertionLength, $change) = @_;
    $$leftFlank = uc($$leftFlank);
    $$rightFlank = uc($$rightFlank);
    $chromSeq =~ m/(.*?$$leftFlank)(.*?)$$rightFlank/ or die "build_chrom.pl error: failed to match change:\n$$change\n";
    my $lastBefore = length($1);   
    my $firstAfter = $lastBefore + length($2) + 1;
    COORDINATES(\$lastBefore, \$firstAfter, $toInsert, $insertionLength);
}

# change lift feature coordinates
sub liftOver {
    my ($lastBefore, $firstAfter, $insertionLength) = @_;
    my $deletionLength = $firstAfter - $lastBefore - 1;     
    my $netInsertionLength = $insertionLength - $deletionLength;
    liftOverPositions($lastBefore, $firstAfter, $insertionLength, $netInsertionLength, \%revLiftStartIs, \@liftStartIs);
    liftOverPositions($lastBefore, $firstAfter, $insertionLength, $netInsertionLength, \%revLiftEndIs, \@liftEndIs);
}
sub liftOverPositions{
    my ($lastBefore, $firstAfter, $insertionLength, $netInsertionLength, $revLifts, $is) = @_;
    foreach my $i(@$is){ 
        my $pos = $$revLifts{$i};
        $pos == -1 and next;
        if($pos >= $firstAfter){  
            $$revLifts{$i} += $netInsertionLength;   
        } elsif($pos > $lastBefore + $insertionLength){  
            $$revLifts{$i} = -1;        
        } else {
            return;
        }
    }
}

# bed file handles
sub getOutBedH {
    my ($bedFile, $type) = @_;
    $type or $type = ">";
    unless($outBedHs{$bedFile}){
        open my $outH, $type, $bedFile or die "build_chrom.pl error: could not open $bedFile for writing: $!\n";
        $outBedHs{$bedFile} = $outH;
    }
    return $outBedHs{$bedFile};
}
sub closeOutBedHs {
    foreach my $bedFile(keys %outBedHs){
        my $outH = $outBedHs{$bedFile};
        close $outH;
    }
    %outBedHs = ();
}

# add new BED features
sub addNewFeatures {
    print "printing requested new $chrom features\n";
    $chromSeq or loadChrom();    
    open my $inH, "-|", $features or die "build_chrom.pl error: could not open features stream for reading: $!\n";
    my $n;
    while(my $feature = <$inH>){
        chomp $feature;  
        $feature =~ s/\r//g;
        $feature or next;
        $feature =~ m/^\s*#/ and next;
        my ($outBedFile, $chrom, $sequence, @trailing) = split(/\s+/, $feature);
        $sequence or next;
        $n++;
        $outBedFile =~ m|^/| or $outBedFile = "$outDir/$outBedFile";    
        $outBedFile .= ".chrom_$chrom";
        my $outH = getOutBedH($outBedFile, ">>");
        $sequence = uc($sequence);
        my $sequenceLength = length($sequence);        
        my $prevEnd = 0;
        my $workingChromSeq = $chromSeq;
        while($workingChromSeq =~ s/(.*?)$sequence//){
            my $start = $prevEnd + length($1);
            my $end = $start + $sequenceLength;
            print $outH join("\t", $chrom, $start, $end, @trailing), "\n";
            $prevEnd = $end;
        }
        $prevEnd or print "build_chrom.pl warning: requested feature line $n had no matches in revised $chrom\n";          
    }
    close $inH;
    closeOutBedHs();
}

