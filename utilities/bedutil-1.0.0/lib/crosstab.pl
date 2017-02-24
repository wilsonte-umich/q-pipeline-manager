#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'crosstab' takes feature BED lines on STDIN as input that have had an
# evaluation score added as the last column, and a simulation iteration
# added as the penultimate column. A crosstable is printed to STDOUT
# in which original input BED6 features are listed in rows, iterations
# are listed in subsequent columns, and the scores are listed at the
# intersection of the associated row/column combinations. The BED name
# field, field 4, is used as the feature unique identifier for making
# row-column associations.  The output format is thus:
#     chrom start end name score strand I-1 I0 I1 ...
#     chr12 123   456 id1  0     +      5   8  3  ...
#     ...
#----------------------------------------------------------------------
# Options are:
#     crosstab has no options
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 20, 14);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
);
%booleanOptions = map { $_ => 1} qw (
);
%defaultValues = (
);
@requiredOptions = qw (
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
$ENV{IS_CHILD} and return 1;
#######################################################################


#######################################################################
# create and print crosstab
#----------------------------------------------------------------------
my (%bed, %features, %iterations, %scores, $isStringID);
while (my $line = <STDIN>){
    my ($strandIndex, $feature) = parseBedLine($line, $error);
    my $maxI = scalar(@$feature) - 1;
    $maxI >= 7 or die "$error: insufficient columns in bed input line:\n$line";
    my $id = $$feature[3];
    $isStringID = $id =~ m/^\d+$/ ? $isStringID : 1;
    my $score = $$feature[$maxI];
    my $iteration = $$feature[$maxI - 1];
    $iteration == -1 and $bed{$id} = join("\t", @$feature[0..5]);
    $features{$id}++;
    $iterations{$iteration}++;
    my $key = "$id,$iteration";
    $scores{$key} and die "$error: encountered more than one score for combination:\nfeature = $id, iteration = $iteration\n";
    $scores{$key} = $score;
}
my @features;
if($isStringID){
    @features = sort {$a cmp $b} keys %features;
} else {
    @features = sort {$a <=> $b} keys %features;
}
my @iterations = sort {$a <=> $b} keys %iterations;
print join("\t", qw(chrom start end name score strand));
foreach my $iteration(@iterations){
    print "\tI$iteration";
}
print "\n";
foreach my $id(@features){
    print $bed{$id};
    foreach my $iteration(@iterations){
        my $key = "$id,$iteration";
        my $score = $scores{$key};
        defined $score or die "$error: missing score for combination:\nfeature = $id, iteration = $iteration\n";
        print "\t$score";
    }
    print "\n";
}
#######################################################################

1;
