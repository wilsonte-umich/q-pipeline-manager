#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'correlate' takes two merged input BED streams and calculates linear 
# correlation coefficients between their score fields.  The features are 
# provided on  STDIN, with the name field (field 4) identifying the source 
# samples. Each sample is expected to have a correspondent set of features, 
# with each feature present at most one time in each sample, e.g.:
#      chr1 0    1000 sample1 55 +
#      chr1 2000 3000 sample1 22 +
#      ...
#      chr1 1000 2000 sample2 33 +
#      chr1 2000 3000 sample2 66 +
#      ...
# If a feature is present in one sample but absent from the other, it is
# assumed to have a score of 0 in the missing sample. The following correlation 
# coefficients are reported to STDOUT, in addition to other summary statistics:
#      base-weighted Pearson (larger features are weighted more heavily)
#      unweighted Pearson (each feature has equal weight regardless of length)
#      unweighted Spearman (based on ranking features by score)
#----------------------------------------------------------------------
# Options are:
#     STRANDED   boolean indicating whether features are strand-specific
#                if stranded, features are compared separately on each strand
#                OPTIONAL [default: FALSE] 
#----------------------------------------------------------------------
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil correlate STRANDED
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 33, 27);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    STRANDED
);
%booleanOptions = map { $_ => 1} qw (
    STRANDED
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
# collect and collapse BED features from STDIN
#----------------------------------------------------------------------
loadBedFile(undef, $error, \my$nFeatures, undef, \my%features);
$nFeatures or exit;
my $names = correlateBED("Pearson", \%features);
setSpearmanRanks(\%features, $names);
correlateBED("Spearman", \%features); 
#######################################################################


#######################################################################
# correlation subs
#----------------------------------------------------------------------
sub correlateBED {
    my ($method, $features) = @_;
    my $isPearson = ($method eq "Pearson");
    my (%names, @names, 
        %wSums,  %wxSums,  %wMeans, 
        %uwSums, %uwxSums, %uwMeans, 
        %wdsSums,  %wStDevs, 
        %uwdsSums, %uwStDevs, 
        $wdpSum,  $wCov,  $wCor,
        $uwdpSum, $uwCov, $uwCor);
#----------------------------------------------------------------------   
    foreach my $chrom(keys %features){  # loop 1: extract and check names
        foreach my $strandIndex(keys %{$features{$chrom}}){
            my $strand = $strandIndex ? ",$strandIndex" : "";
            foreach my $start(keys %{$features{$chrom}{$strandIndex}}){   
                foreach my $end(keys %{$features{$chrom}{$strandIndex}{$start}}){
                    my $trailings = $features{$chrom}{$strandIndex}{$start}{$end};
                    my $nTrailings = @$trailings;
                    $nTrailings > 2 and die "$error: $nTrailings features present for $chrom:$start-$end$strand: expected 1 or 2\n";                     
                    foreach my $trailing(@$trailings){
                        my ($name) = @$trailing;  
                        $names{$name}++;
                    }  
                }
            }     
        }
    }    
    @names = sort {$a cmp $b} keys %names;
    my $nNames = @names;
    $nNames == 2 or die "$error: $nNames samples present in input stream: expected 2\n";                
#----------------------------------------------------------------------
    foreach my $chrom(keys %features){  # loop 2: collect scores
        foreach my $strandIndex(keys %{$features{$chrom}}){
            foreach my $start(keys %{$features{$chrom}{$strandIndex}}){   
                foreach my $end(keys %{$features{$chrom}{$strandIndex}{$start}}){
                    my $weight = $end - $start;
                    my %encountered;
                    foreach my $trailing(@{$features{$chrom}{$strandIndex}{$start}{$end}}){
                        my ($name, $score, $strand) = @$trailing;  
                        $wSums{$name} += $weight;
                        $wxSums{$name} += ($weight * $score);
                        $uwSums{$name}++;
                        $uwxSums{$name} += $score;
                        $encountered{$name}++;
                    }  
                    foreach my $name(@names){
                        unless($encountered{$name}){
                            $wSums{$name} += $weight;
                            $uwSums{$name}++;
                        }
                    }   
                }
            }     
        }
    }       
#----------------------------------------------------------------------    
    foreach my $name(@names){   # check data integrity and calculate mean
        $wSums{$name} or die "$error: summed feature length for sample $name was zero\n";
        $wMeans{$name} = $wxSums{$name} / $wSums{$name};
        $uwMeans{$name} = $uwxSums{$name} / $uwSums{$name};
    }         
#----------------------------------------------------------------------    
    foreach my $chrom(keys %features){  # loop 3: collect delta values dependent on mean
        foreach my $strandIndex(keys %{$features{$chrom}}){
            foreach my $start(keys %{$features{$chrom}{$strandIndex}}){   
                foreach my $end(keys %{$features{$chrom}{$strandIndex}{$start}}){
                    my $weight = $end - $start;
                    my (%wDeltas, %uwDeltas);
                    my %encountered;
                    foreach my $trailing(@{$features{$chrom}{$strandIndex}{$start}{$end}}){
                        my ($name, $score, $strand) = @$trailing;  
                        $wDeltas{$name} = $score - $wMeans{$name};
                        $uwDeltas{$name} = $score - $uwMeans{$name};
                        $wdsSums{$name} += $weight * $wDeltas{$name} ** 2;
                        $uwdsSums{$name} += $uwDeltas{$name} ** 2;    
                        $encountered{$name}++;              
                    }
                    foreach my $name(@names){
                        unless($encountered{$name}){
                            $wDeltas{$name} = -$wMeans{$name};
                            $uwDeltas{$name} = -$uwMeans{$name};
                        }
                    } 
                    $wdpSum += $weight * $wDeltas{$names[0]} * $wDeltas{$names[1]};
                    $uwdpSum += $uwDeltas{$names[0]} * $uwDeltas{$names[1]};
                }
            }  
        }
    } 
#----------------------------------------------------------------------  
    foreach my $name(@names){  # calculate stdev, covariance and correlation coefficient
        $wStDevs{$name} = sqrt( $wdsSums{$name} / ($wSums{$name} - 1) );
        $uwStDevs{$name} = sqrt( $uwdsSums{$name} / ($uwSums{$name} - 1) );  
    }  
    $wCov = $wdpSum / ($wSums{$names[0]} - 1);
    $wCor = $wCov / ( $wStDevs{$names[0]} * $wStDevs{$names[1]} );
    $uwCov = $uwdpSum / ($uwSums{$names[0]} - 1);
    $uwCor = $uwCov / ( $uwStDevs{$names[0]} * $uwStDevs{$names[1]} );
#----------------------------------------------------------------------  
    my $divider = "-" x 50;  # print the results
    if($isPearson){
        print "$divider\nGeneral\n$divider\n";
        print "$uwSums{$names[0]}\tfeatures in each feature set\n";
        print "$wSums{$names[0]}\tbases in each feature set\n";
        foreach my $name(@names){
            print "$uwMeans{$name}\t$name\tunweighted mean\n";
            print "$uwStDevs{$name}\t$name\tunweighted stddev\n";        
            print "$wMeans{$name}\t$name\tweighted mean\n";
            print "$wStDevs{$name}\t$name\tweighted stddev\n";    
        }         
    }
    print "$divider\n$method\n$divider\n";    
    print "$uwCov\tunweighted covariance\n";
    print "$uwCor\tunweighted $method correlation coefficient\n";  
    $isPearson or (print "$divider\n" and return);  # weighted Spearman is not valid as bases were never ranked, only features
    print "$wCov\tweighted covariance\n";
    print "$wCor\tweighted $method correlation coefficient\n"; 
#---------------------------------------------------------------------- 
    return \@names;   
} 
#----------------------------------------------------------------------
sub setSpearmanRanks {  # convert raw scores to Spearman ranks
    my ($features, $names) = @_;
    my %scores;
#----------------------------------------------------------------------
    foreach my $chrom(keys %features){  # loop 4: collect scores for ranking
        foreach my $strandIndex(keys %{$features{$chrom}}){
            foreach my $start(keys %{$features{$chrom}{$strandIndex}}){   
                foreach my $end(keys %{$features{$chrom}{$strandIndex}{$start}}){
                    my %encountered;
                    foreach my $trailing(@{$features{$chrom}{$strandIndex}{$start}{$end}}){
                        my ($name, $score, $strand) = @$trailing; 
                        push @{$scores{$name}{$score}}, [$chrom, $strandIndex, $start, $end];
                        $encountered{$name}++;
                    } 
                    foreach my $name(@$names){
                        unless($encountered{$name}){
                            push @{$scores{$name}{0}}, [$chrom, $strandIndex, $start, $end];
                        }
                    }   
                }
            }  
        }
    } 
#----------------------------------------------------------------------
    my $i = 1;
    %$features = ();
    foreach my $name(keys %scores){  # convert scores to Spearman ranks
        foreach my $score(sort {$a <=> $b} keys %{$scores{$name}}){
            my $scores = $scores{$name}{$score};
            my $nScores = @$scores;
            my $rank = getSpearmanRank($i, $nScores);
            foreach my $feature(@$scores){
                my ($chrom, $strandIndex, $start, $end) = @$feature;
                push @{$$features{$chrom}{$strandIndex}{$start}{$end}}, [$name, $rank, $strandIndex?$strandIndex:"+"]
            }
            $i += $nScores;
        }
    } 
}
#----------------------------------------------------------------------
sub getSpearmanRank {
    my ($i, $nScores) = @_;
    $nScores == 1 and return $i;
    my @ranks = ($i..($i + $nScores - 1));
    my $rankSum;
    $rankSum += $_ for @ranks;
    return $rankSum / $nScores;
}
#######################################################################

1;

