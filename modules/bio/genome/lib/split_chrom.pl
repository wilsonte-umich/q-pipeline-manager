#!/usr/bin/perl
use strict;
use warnings;
use Cwd(qw(abs_path));

#q    require $GENOME_DIR $GENOME $IS_FASTA $OUT_FILE $CHROMOSOME $SPLIT_SEQ $NAME

#$    -N  split.$GENOME.$CHROMOSOME.$NAME
#$    -wd $GENOME_DIR
#$    -l  vf=4G
#$    -l  h_rt=5:00:00

#PBS  -N  split.$GENOME.$CHROMOSOME.$NAME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=4gb
#PBS  -l  walltime=5:00:00

# define variables
my $genomeDir = $ENV{GENOME_DIR};
my $genome = $ENV{GENOME};
my $chrom = $ENV{CHROMOSOME};
my $splitSeq = uc($ENV{SPLIT_SEQ});
my $rcSplitSeq = rc($splitSeq);
$splitSeq   =~ s/N/\./g; # allow wildcard matches
$rcSplitSeq =~ s/N/\./g;
my $seqName = $ENV{NAME};
print "splitting $genome $chrom on sequence $seqName\n";
my $chromFasta = "$genomeDir/chromosomes/$chrom.fa";
-e $chromFasta or die "split_chrom.pl error: file not found: $chromFasta\n";
my $isFasta = $ENV{IS_FASTA};
my $outFile = "$ENV{OUT_FILE}.chrom_$chrom";
my ($prevEnd, $chromSeq, $before, $seq, $outH) = 0;

# main execution block
loadChrom();
splitChrom();
print "done\n";

# worker subs
sub rc {
    my ($seq) = @_;
    $seq =~ tr/ACGT/TGCA/;
    return reverse($seq);
}
sub loadChrom {
    open my $inH, "<", $chromFasta or die "split_chrom.pl error: could not open $chromFasta for reading: $!\n";
    while(my $line = <$inH>){
        chomp $line;
        $line =~ s/\s//g;
        $line or next;
        $line =~ m/^>/ and next;
        $chromSeq .= uc($line);
    }
    close $inH;
}
sub splitChrom {
    open $outH, ">", $outFile or die "split_chrom.pl error: could not open $outFile for writing: $!\n";
    while($chromSeq =~ s/(.*?)($splitSeq|$rcSplitSeq)//){
        $before = $1;
        $seq = $2;
        $prevEnd = commitSeq(); 
    }
    $before = $chromSeq;
    $seq = "";
    commitSeq();
    close $outH;
}
sub commitSeq {
    ($before or $seq) or return $prevEnd;
    my $start = $prevEnd;
    my $length = length($before) + length($seq);
    my $end = $start + $length;
    my $start1 = $start + 1;
    if ($isFasta) {
        my $name = "$chrom\_$start1\_$end";
        print $outH ">$name\n$before$seq\n";
    } else {
        my $name = "$chrom:$start1-$end";            
        my $strand = $seq eq $splitSeq ? "+" : "-";
        print $outH join("\t", $chrom, $start, $end, $name, $length, $strand), "\n";  
    }
    return $end;
}

