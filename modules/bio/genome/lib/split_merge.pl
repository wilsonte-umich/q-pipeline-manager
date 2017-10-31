#!/usr/bin/perl
use strict;
use warnings;

use Cwd(qw(abs_path));
#q    require $GENOME_DIR $GENOME $OUT_FILE $NAME $IS_FASTA

#$    -N  split.$GENOME.merge.$NAME
#$    -wd $GENOME_DIR
#$    -l  vf=2G

#PBS  -N  split.$GENOME.merge.$NAME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=2gb

# define variables
my $genomeDir = $ENV{GENOME_DIR};
my $genome = $ENV{GENOME};
my $seqName = $ENV{NAME};
print "merging $genome splits on sequence $seqName\n";
my $outFile = $ENV{OUT_FILE};
my $inGlob = "$outFile.chrom_*";
my $isFasta = $ENV{IS_FASTA};
my ($minLength, $maxLength) = (1e9, 0);
my (%lengths, %bins);

# merge
print "$outFile\n";
qx/slurp -o $outFile $inGlob/;
qx/rm -f $inGlob/;
print qx/head $outFile; echo/;

# collect lengths
open my $inH, "<", $outFile or die "could not open $outFile for reading: $!\n";
while (my $line = <$inH>) {
    chomp $line;
    $line or next;
    my $length;
    if ($isFasta) {
        $line =~ m/^>/ or next;
        my @bits = split("_", $line);
        my $end = pop @bits;
        my $start = pop @bits;
        $length = $end - $start + 1;
    } else {
        my @bits = split("\t", $line);
        $length = $bits[4];
    }
    $lengths{$length}++;
    $minLength <= $length or $minLength = $length;
    $maxLength >= $length or $maxLength = $length;
}
close $inH;

# print histogram
my $spread = $maxLength - $minLength;
my $binSize = int($spread / 100 + 0.5);
foreach my $length (keys %lengths) {
    my $bin = int($length / $binSize + 0.5) * $binSize;
    $bins{$bin} += $lengths{$length};
}
print "length_bin\tcount\n";
foreach my $bin (sort {$a <=> $b} keys %bins) {
    print "$bin\t$bins{$bin}\n";
}

print "done\n";

