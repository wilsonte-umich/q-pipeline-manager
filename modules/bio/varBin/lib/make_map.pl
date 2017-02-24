use strict;
use warnings;

# takes a list of pos, [lens] and converts it to a set of coverage break points

my ($CHROM_SIZE, $STAT_FILE) = @ARGV;
my ($N, $sLen, $pos0, @breaks) = (0, 0, 0);

# parse the break points
while(<STDIN>){
    chomp;
    my ($pos, $lens) = split("\t");
    my @lens = split(",", $lens);
    $pos0 and printBreaks($#breaks < $pos-$pos0-1 ? $#breaks : $pos-$pos0-1);
    $pos0 = $pos; 
    $breaks[0] += @lens; # increment the coverage plot up
    foreach my $len(@lens){
        $sLen += $len;
        $N++;
        $breaks[$len]--; # increment the coverage plot down
    }
}
printBreaks($#breaks);

open(my $statH, ">", $STAT_FILE) or die " could not open $STAT_FILE for writing: $!\n";
print $statH join("\t", 'nFrags',   $N), "\n";
print $statH join("\t", 'sumLen',   $sLen), "\n";
print $statH join("\t", 'chromLen', $CHROM_SIZE), "\n";
print $statH join("\t", 'avgLen',   int($sLen / $N + 0.5)), "\n";
print $statH join("\t", 'coverage', int($sLen / $CHROM_SIZE * 10 + 0.5) / 10), "\n";
close $statH;

sub printBreaks {
    my ($maxI) = @_;
    foreach my $i(0..$maxI){
        my $n = shift @breaks;
        $n or next;
        print pack('L', $i + $pos0). # position as unsigned 4 byte/32-bit = 4,294,967,295
              pack('s', $n);         # increment value as signed 2 byte/16-bit = 32,767           
    }  
}
