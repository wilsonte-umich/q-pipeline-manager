#!/usr/bin/perl
use strict;
use warnings;
use Cwd(qw(abs_path));

#q    require $INPUT_GENOME_DIR $INPUT_GENOME $LIFT_BEDS $OUTPUT_GENOME_DIR $OUTPUT_GENOME $LIFT_BEDS $FEATURES_FILE

#$    -N  build.$OUTPUT_GENOME.merge
#$    -wd $GENOMES_DIR
#$    -l  vf=2G
#$    -l  h_rt=5:00:00

#PBS  -N  build.$OUTPUT_GENOME.merge
#PBS  -d  $GENOMES_DIR
#PBS  -l  mem=2gb
#PBS  -l  walltime=5:00:00

# define variables
my $inGenome = $ENV{INPUT_GENOME};
my $outGenome = $ENV{OUTPUT_GENOME};
print "merging $outGenome chromosome files\n\n";
my $inDir = $ENV{INPUT_GENOME_DIR};
my $outDir = $ENV{OUTPUT_GENOME_DIR};
my $chromsFasta = "$outDir/chromosomes/*.fa";
my $outFasta = "$outDir/$outGenome.fa";
my @inLiftBeds = getInLiftBeds($ENV{LIFT_BEDS});
my @featureBeds = getFeatureBeds($ENV{FEATURES_FILE});
my %bedFiles = map { $_ => 1 } (@inLiftBeds, @featureBeds);
my $statsFile = "$outDir/$outGenome.stats.csv";

# main execution block
print "$outFasta\n";
qx/slurp -o $outFasta $chromsFasta/;
print qx/head $outFasta; echo/;
print qx/samtools faidx $outFasta/;
my $nBases = "slurp $outFasta | grep -vP '^>' | sed 's/\\s//g' | awk '{c+=length(\$0)}END{print c}'";
$nBases = qx/$nBases/;
chomp $nBases;
qx/set_stat $outDir $outGenome N_BASES,$nBases/;
foreach my $bedFile(keys %bedFiles){
    print "$bedFile\n";
    my $inGlob = "$bedFile.chrom_*";
    qx/slurp -o $bedFile $inGlob/;
    qx/rm -f $inGlob/;
    print qx/head $bedFile; echo/;
}
print "$statsFile\n";
print qx/cat $statsFile; echo/;
print "done\n";

# bed file subs
sub getInLiftBeds {
    my ($inLeftBeds) = @_;
    $inLeftBeds or return ();
    return map { getOutLiftBed($_) } split(" ", $inLeftBeds); 
}
sub getOutLiftBed {
    my ($inLiftBed) = @_;
    -e $inLiftBed or $inLiftBed = "$inDir/$inLiftBed";
    my $outLiftBed = $inLiftBed;
    $outLiftBed =~ s/$inGenome/$outGenome/g;
    return $outLiftBed;
}
sub getFeatureBeds {
    my ($featuresFile) = @_;
    $featuresFile or return ();
    open my $inH, "<", $featuresFile or die "build_merge.pl error: could not open $featuresFile for reading: $!\n";
    my %outBedFiles;
    while(my $feature = <$inH>){
        chomp $feature;  
        $feature =~ s/\r//g;
        $feature or next;
        $feature =~ m/^\s*#/ and next;
        my ($outBedFile) = split(/\s+/, $feature); 
        $outBedFile =~ m|^/| or $outBedFile = "$outDir/$outBedFile";      
        $outBedFiles{$outBedFile}++;        
    }
    close $inH; 
    return keys %outBedFiles;
}

