#!/usr/bin/perl

#q    require $GENOME_DIR $GENOME $ANNOT $TRANS_TXT 
#q    require $NAMES_TXT $TRANS_ID_COL $GENE_ID_COL
#q    require $MAP_STR_BED $MAP_UNSTR_BED $MAP_GENES_BED 

#$    -N  map_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=4G
#$    -l  h_rt=2:00:00

#PBS  -N  map_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=4gb
#PBS  -l  walltime=2:00:00

#===============================================================================
#synopsis
#-------------------------------------------------------------------------------
#map.pl takes UCSC-formatted transcript isoforms, and creates a BED map of the transcriptome
#such that each and every genome position (and sometimes strand) is assigned into a
#single segment of one of the following types:
#   exon                (annotated exon entries in unstranded file)   
#   exon_sense          (annotated exon entries in stranded file)
#   exon_antisense      (reflection of exon entries onto the opposite strand)
#   intron              (annotated intron entries in unstranded file)   
#   intron_sense        (annotated intron entries in stranded file)
#   intron_antisense    (reflection of intron entries onto the opposite strand)
#   ambiguous           (regions that cannot be resolved to a single preferred element)  
#   intergenic          (regions between annotated transcription elements; not included)
#Two output BED files are created: an 'unstranded' file that contains 'exon' and
#'intron' entries for which overlap is determined without respect to strand, and a
#'stranded' file that contains 'exon_sense' and similar entries for which overlap
#is determined with respect to strand, e.g. exons on different strands are not
#considered to be overlapping.  Choice of the identity at any given genome position 
#(and sometimes strand) is determined by preferring, in order:
#     sense elements over anti-sense elements (stranded file only)
#     exonic elements over intronic elements
#During this process, transcript isoforms of the same gene are condensed into a
#singular "gene".  When an exon is assigned to more than one gene, a purging process
#eliminates common known patterns of duplicated entries where different names
#refer to the same gene.  When these rules do not lead to a singular assignment 
#at a series of contiguous positions (on a strand), the segment is called 'ambiguous'.
#Finally, a BED file of the assembled gene spans in created.
#===============================================================================

#===============================================================================
#OUTPUT FILE format, extension of BED
#-------------------------------------------------------------------------------
#  1  chrom
#  2  start (0-based)
#  3  end (1-based, i.e. half-open)
#  4  name (gene identifier or na)
#  5  score (per BED, 0 to 1000, only used for display purposes)
#  6  strand (+ or -)
#-------------------------------------------------------------------------------
#  7  gene length (span of gene, including all exons)
#  8  feature type (e.g. exon, intron_antisense, gene)
#  9  splice index in gene (left-to-right order, 1-based, 0=na)
#===============================================================================

#===============================================================================
#main program block
#-------------------------------------------------------------------------------
use strict;
use warnings;
my $i = 0;
my @fields = $ENV{ANNOT} eq "UCSC" ? qw (
    name
    chrom
    strand
    txStart
    txEnd
    cdsStart
    cdsEnd
    exonCount
    exonStarts
    exonEnds
    proteinID
    alignID  
) : qw (
    bin
    name
    chrom
    strand
    txStart
    txEnd
    cdsStart
    cdsEnd
    exonCount
    exonStarts
    exonEnds
    score
    name2
    cdsStartStat
    cdsEndStat
    exonFrames
);
my %fields = map { $_ => $i++ } @fields;  	
defined $fields{name2} or $fields{name2} = $fields{name};  # UCSC data lacks name2, use name, but should be overruled by geneSymbol lookup
my $transcriptsFile = $ENV{TRANS_TXT};  # transcripts file as supplied by UCSC
-e $transcriptsFile or (print "$transcriptsFile does not exist\nexiting quietly\n" and exit);
my $geneSymbolsFile = $ENV{NAMES_TXT};  # file for converting transcript IDs to gene symbols
-e $geneSymbolsFile or (print "$geneSymbolsFile does not exist\nexiting quietly\n" and exit);
my (%geneSymbols, %transcripts, %nameSpans, %geneSpans, %exonSpans, %revNames, %exons, 
    %blocks, %geneLengths, $overlapBlock);
my $unStrn = '.';
my $delimiter = "~~~";
$| = 1;
loadGeneSymbols();
getTranscripts();
findNameSpans();
%transcripts = ();
parseExons();
%nameSpans = ();
printGenesBed();
%geneSpans = ();
purgeUnwantedDups();
%exonSpans = ();
parseIntrons();
parseSegments();
print "done\n";
#===============================================================================

#===============================================================================
sub loadGeneSymbols { #load conversion of transcriptIDs to geneSymbols
    print "loading $ENV{GENOME} $ENV{ANNOT} transcriptID to geneSymbol conversions\n";
    open my $inH, "<", $geneSymbolsFile or die "could not open $geneSymbolsFile for reading: $!\n";
    my $n = 0;
    while (my $line = <$inH>){
        chomp $line;
        my @fields = split("\t", $line);
        my $transcriptID = $fields[$ENV{TRANS_ID_COL}-1];
        my $geneSymbol =   $fields[$ENV{GENE_ID_COL}-1];        
        $transcriptID =~ m/\w/ or next;
        $geneSymbol =~ m/\w/ or next;        
        $geneSymbols{$transcriptID} = $geneSymbol;
        $n++;
    }       
    close $inH;  
    print "  $n transcript IDs loaded\n";
}
#-------------------------------------------------------------------------------
#get the unique gene spans and exons for each chrom/strand/gene combination
#discard splice variant information
sub getTranscripts { # load the transcripts from the transcripts file
    print "loading $ENV{GENOME} $ENV{ANNOT} transcripts\n";
    open my $inH, "<", $transcriptsFile or die "could not open $transcriptsFile for reading: $!\n";
    my $n = 0;
    while (my $line = <$inH>){
        chomp $line;
        my @fields = split("\t", $line);
        my $chrom = $fields[$fields{chrom}];
        $chrom =~ m/_/ and next; #only map the canonical chromosomes    
        my $transcriptID = $fields[$fields{name}];   
        $transcriptID or next;
        my $name = $geneSymbols{$transcriptID} ? $geneSymbols{$transcriptID} : $fields[$fields{name2}];
        $name or next;  # from here forward, name is the geneSymbol, from lookup file or name2 field, in that precedence order
        #$name =~ m/\-AS\d+$/ and next; #ignore annotated antisense transcripts
        $n++;
        my $strand = $fields[$fields{strand}];
        my @exonStarts = split(",", $fields[$fields{exonStarts}]); 
        @exonStarts = map { $_ + 1 } @exonStarts; # use 1-referenced starts in this script; convert back on printing
        my @exonEnds = split(",", $fields[$fields{exonEnds}]);    
        my $txnStart = $fields[$fields{txStart}] + 1;
        my $txnEnd = $fields[$fields{txEnd}];
        push @{$transcripts{$name}{$chrom}{$strand}}, [$txnStart, $txnEnd, \@exonStarts, \@exonEnds]; 
    }
    close $inH;  
    print "  $n transcripts loaded from canonical chromosomes\n";
}
#-------------------------------------------------------------------------------
sub findNameSpans {  # find the non-overlapping genome spans covered by transcripts of each named gene
    print "finding genome spans for named genes\n";
    my ($nNames, $nNameSpans);
    foreach my $name(keys %transcripts){
        $nNames++;
        foreach my $chrom(keys %{$transcripts{$name}}){
            foreach my $strand(keys %{$transcripts{$name}{$chrom}}){
                my $transcripts = $transcripts{$name}{$chrom}{$strand};
                my $nTranscripts = @$transcripts;
                if($nTranscripts == 1){  # only one span on this chromosome strand, no need to do overlap discovery
                    my ($txnStart, $txnEnd, $exonStarts, $exonEnds) = @{$$transcripts[0]};
                    my $spanIndex = join($delimiter, $chrom, $strand, $txnStart, $txnEnd);
                    @{$nameSpans{$name}{$spanIndex}} = ([$exonStarts, $exonEnds]);  
                    $nNameSpans++;
                } else { 
                    my $overlapGroups = findOverlapGroups($transcripts);
                    $nNameSpans += commitOverlapGroups($overlapGroups, $name, $chrom, $strand);
                }    
            } 
        }
    }
    print "  $nNames gene names encountered, mapping to $nNameSpans unique genome spans\n";
}
sub findOverlapGroups {  # discover whether multiple transcripts are overlapping, if not, split into separate spans
    my ($transcripts) = @_;
    my %overlapGroups;
    TRANSCRIPT: foreach my $transcript(@$transcripts){
        my ($txnStart, $txnEnd, $exonStarts, $exonEnds) = @$transcript;
        my $outTranscripts = [[$exonStarts, $exonEnds]];
        collapseTranscriptsSpan(\%overlapGroups, $outTranscripts, $txnStart, $txnEnd); 
    } 
    scalar(keys %overlapGroups) > 1 and collapseOverlapGroups(\%overlapGroups);
    return \%overlapGroups;
}
sub collapseTranscriptsSpan {  # process a span against any current overlap groups; collapse it once, if possible
    my ($overlapGroups, $inTranscripts, $inStart, $inEnd) = @_;
    if(scalar(keys %$overlapGroups)){
        foreach my $overlapIndex(keys %$overlapGroups){
            my $groupTranscripts = $$overlapGroups{$overlapIndex};
            my ($groupStart, $groupEnd) = split($delimiter, $overlapIndex);
            if($groupStart <= $inEnd and $inStart <= $groupEnd){  # overlaps previous group, add to it
                delete $$overlapGroups{$overlapIndex}; # reset group boundaries as needed
                $groupStart <= $inStart or $groupStart = $inStart;
                $groupEnd >= $inEnd or $groupEnd = $inEnd;
                @{$$overlapGroups{"$groupStart$delimiter$groupEnd"}} = (@$groupTranscripts, @$inTranscripts);
                return;
            }   
        }
        @{$$overlapGroups{"$inStart$delimiter$inEnd"}} = (@$inTranscripts); #no overlap, is a new group
    } else {  # first encountered span defines the first overlap group
        @{$$overlapGroups{"$inStart$delimiter$inEnd"}} = (@$inTranscripts);
    }
}
sub collapseOverlapGroups { # recursively collapse overlap groups until no more group merging occurs
    my ($overlapGroups) = @_;
    my %newOverlapGroups;
    foreach my $overlapIndex(keys %$overlapGroups){
        my $groupTranscripts = $$overlapGroups{$overlapIndex};
        my ($groupStart, $groupEnd) = split($delimiter, $overlapIndex);
        collapseTranscriptsSpan(\%newOverlapGroups, $groupTranscripts, $groupStart, $groupEnd);
    }
    scalar(keys %$overlapGroups) == scalar(keys %newOverlapGroups) or collapseOverlapGroups(\%newOverlapGroups);
    %$overlapGroups = %newOverlapGroups;
}
sub commitOverlapGroups {  # commit transcript overlap groups as spans corresponding to named gene; might be more than one!
    my ($overlapGroups, $name, $chrom, $strand) = @_;
    my $n;
    foreach my $overlapIndex(keys %$overlapGroups){
        my $groupTranscripts = $$overlapGroups{$overlapIndex};
        my ($groupStart, $groupEnd) = split($delimiter, $overlapIndex);
        my $spanIndex = join($delimiter, $chrom, $strand, $groupStart, $groupEnd);
        @{$nameSpans{$name}{$spanIndex}} = @$groupTranscripts;
        $n++;
    }
    return $n;
}
#-------------------------------------------------------------------------------
sub parseExons {  # parse the exons for each genome span for each named gene
    print "parsing exons out of gene spans\n";
    foreach my $name(keys %nameSpans){
        my $spans = $nameSpans{$name};
        my $nSpans = scalar(keys %$spans);
        my $i = 1;        
        foreach my $spanIndex(sort {$a cmp $b} keys %$spans){
            my $transcripts = $nameSpans{$name}{$spanIndex};
            my ($chrom, $strand, $spanStart, $spanEnd) = split($delimiter, $spanIndex);
            my $splitName = $nSpans == 1 ? $name : "$name\:$i";  # modify the gene name when it corresponds to more than one genome span
            $revNames{$splitName}++;  # remember gene names prior to further processing steps below  
            $geneSpans{$spanIndex}{$splitName}++; # reverse to now record the name(s) assigned to each unique span
            foreach my $transcript(@$transcripts){
                my ($exonStarts, $exonEnds) = @$transcript;
                my $nExons = @$exonStarts;
                foreach my $j(1..$nExons){
                    my $k = $j - 1;
                    my $start = $$exonStarts[$k];
                    my $end = $$exonEnds[$k];
                    my $exonIndex = join($delimiter, $chrom, $start, $end, $strand);
                    $exonSpans{$exonIndex}{$splitName}++;  # record the gene names that map to each unique exon span
                }
            }
            $i++;
        }
    }
    my $nExonSpans = scalar(keys %exonSpans);
    print "  $nExonSpans unique exon spans encountered\n";
}
#===============================================================================

#===============================================================================
#Use a rule-based approach to eliminate common known duplicated gene names
#mapping to a single span/exon.  Typically, these arise not from biology but 
#from a gene being annotated with both a meaningful name as well as a generic 
#name (e.g. 'LOC####') or when two gene family members map to the same 
#coordinates (e.g. XYZ1 and XYZ2).
#-------------------------------------------------------------------------------
sub printGenesBed {
    print "purging unwanted duplicate gene names and printing genes bed\n";
    open my $genesH, ">", $ENV{MAP_GENES_BED} or die "could not open $ENV{MAP_GENES_BED} for writing: $!\n";
    foreach my $spanIndex(keys %geneSpans){
        my ($chrom, $strand, $spanStart, $spanEnd) = split($delimiter, $spanIndex);    
        my @names = keys %{$geneSpans{$spanIndex}};
        @names = purgeDupNames(@names);
        my $name = join(",", @names);  # for genes, concatenate whatever names remain; may result in same name for different spans
        print $genesH join("\t", $chrom, $spanStart - 1, $spanEnd, $name, 0, $strand), "\n";
        #record gene lengths from parsed spans (not by recalculating from most widely separated exons)
        #name may not necessarily correspond to the name arrived at separately for each exon, below
        $geneLengths{$chrom}{$name} = $spanEnd - $spanStart + 1;    
    }
    close $genesH;
}
sub purgeUnwantedDups {
    print "purging unwanted duplicate exon names\n";
    foreach my $span(keys %exonSpans){
        my ($chrom, $start, $end, $strand) = split($delimiter, $span);    
        my @names = keys %{$exonSpans{$span}};
        @names = purgeDupNames(@names);
        foreach my $name(@names){   
            $exons{$chrom}{$strand}{$name}{$start}{$end} = [$strand, $name];
            $exons{$chrom}{$unStrn}{$name}{$start}{$end} = [$strand, $name];
        }
    }
}
#-------------------------------------------------------------------------------
sub purgeDupNames {  # apply the name dup purging rules
    my (@names) = @_;
    if(scalar(@names) > 1){
        @names = stripCompoundNames(@names);
        scalar(@names) > 1 and @names = stripGenericNames('^LOC\d+\:\d+$', @names); # :\d was added if split into groups
        scalar(@names) > 1 and @names = stripGenericNames('^LOC\d+$', @names);
        scalar(@names) > 1 and @names = stripGenericNames('^C.{1,2}orf\d+\:\d+$', @names);
        scalar(@names) > 1 and @names = stripGenericNames('^C.{1,2}orf\d+$', @names);
        scalar(@names) > 1 and @names = stripGenericNames('^Gm\d+\:\d+$', @names);   
        scalar(@names) > 1 and @names = stripGenericNames('^Gm\d+$', @names);    
        scalar(@names) > 1 and @names = condenseGeneGroups(@names);  
        scalar(@names) or die "map.pl error: purgeDupNames resulted in no names for: @names\n";
    }  
    return @names;
}
sub stripCompoundNames { #remove compound gene names (fusion transcripts) where one portion of name is already known as a gene
    my (@names) = @_;
    my @strippedNames;
    foreach my $name(@names){
        ($name =~ m/^(\w+)\-(\w+)$/ and ($revNames{$1} or $revNames{$2})) or push @strippedNames, $name;
    }
    scalar(@strippedNames) and return @strippedNames;
    return @names; # all names were compound, cannot purge
}
sub stripGenericNames { #remove gene names that are of a generic type when at least one non-generic name remains
    my ($regex, @names) = @_;
    my @strippedNames;
    foreach my $name(@names){ $name =~ m/$regex/ or push @strippedNames, $name }
    scalar(@strippedNames) and return @strippedNames;
    return @names; # all names were of the generic type, cannot purge
}
sub condenseGeneGroups { #condense gene groups such as AMY1A,AMY1B or MIR19-1,MIR19-2 into their common root, i.e. AMY1 or MIR19
    my (@names) = @_;
    @names = sort {length($a) <=> length($b)} @names;
    checkGroupRoot($names[0], '', '\w', '+', @names[1..$#names]) and return $names[0];
    my $root;    
    checkGroupForRoot(\$root, '', '[A-Z]', '{1}', @names) and return $root;    
    checkGroupForRoot(\$root, '-', '\d', '{1,2}', @names) and return $root; 
    checkGroupForRoot(\$root, '', '\d', '{1,2}', @names) and return $root;
    my $length = length($names[0]);
    my $minLength = $length - 4;
    $minLength = $minLength > 2 ? $minLength : 2;
    my @lengths = reverse($minLength..$length);
    for my $length(@lengths) { 
        my $root = substr($names[0], 0, $length);
        stepBackForRoot($root, @names[1..$#names]) and return $root;  
    }
    return @names; 
}
sub checkGroupRoot {
    my ($root, $linker, $suffix, $nIter, @names) = @_;
    foreach my $name(@names){ $name =~ m/^$root$linker($suffix)$nIter$/ or return 0 }
    return 1; 
}
sub checkGroupForRoot {
    my ($root, $linker, $suffix, $nIter, @names) = @_;
    if($names[0] =~ m/^(.+?)$linker($suffix)$nIter$/){ 
        $$root = $1;
        return checkGroupRoot($$root, $linker, $suffix, $nIter, @names[1..$#names]);
    }
}
sub stepBackForRoot {
    my ($root, @names) = @_;
    foreach my $name(@names){ $name =~ m/^$root/ or return 0 }
    return 1;   
}
#===============================================================================

#===============================================================================
#merge to the widest possible exons, and the introns between them, WITHIN each chrom/strand/gene combination
#-------------------------------------------------------------------------------
sub parseIntrons { 
    print "parsing introns\n";
    foreach my $chrom(sort {$a cmp $b} keys %exons){
        foreach my $strandIndex(sort {$a cmp $b} keys %{$exons{$chrom}}){
            foreach my $gene(sort {$a cmp $b} keys %{$exons{$chrom}{$strandIndex}}){
                my ($exonCounter, $intronCounter, $exonStart, $exonEnd, $fields, $prevExonEnd) = (0, 0, 0, 0); 
                foreach my $start(sort {$a <=> $b} keys %{$exons{$chrom}{$strandIndex}{$gene}}){
                    if($exonEnd){
                        if($start > $exonEnd){
                            commitSplice(\$exonCounter, \$intronCounter, $chrom, $strandIndex, $exonStart, $exonEnd, $fields, $prevExonEnd);
                            $prevExonEnd = $exonEnd;
                            $exonStart = $start;
                            $exonEnd = 0;
                        }  
                    } else {
                        $exonStart = $start;
                    }  
                    foreach my $end(sort {$a <=> $b} keys %{$exons{$chrom}{$strandIndex}{$gene}{$start}}){                
                        $exonEnd = $exonEnd >= $end ? $exonEnd : $end;
                        $fields = $exons{$chrom}{$strandIndex}{$gene}{$start}{$end};
                    }
                }    
                commitSplice(\$exonCounter, \$intronCounter, $chrom, $strandIndex, $exonStart, $exonEnd, $fields, $prevExonEnd);
            }
        }
    }
}
sub commitSplice {
    my ($exonCounter, $intronCounter, $chrom, $strandIndex, $start, $end, $fields, $prevEnd) = @_;
    my ($strand, $name) = @$fields;
    $$exonCounter = $$exonCounter ? $$exonCounter += 2 : 1;
    my $spliceIndex = $$exonCounter;
    if($strandIndex eq $unStrn){
        push @{$blocks{$chrom}{$strandIndex}{$start}{$end}}, ['exon', $strand, $name, $spliceIndex];
    } else {
        my $antisenseStrandIndex = $strandIndex eq '+' ? '-' : '+';
        push @{$blocks{$chrom}{$strandIndex}         {$start}{$end}}, ['exon_sense', $strandIndex, $name, $spliceIndex];
        push @{$blocks{$chrom}{$antisenseStrandIndex}{$start}{$end}}, ['exon_antisense', $antisenseStrandIndex, $name, $spliceIndex];
    }
    if($prevEnd){
        $$intronCounter = $$exonCounter - 1;
        $spliceIndex = $$intronCounter;
        my $intronStart = $prevEnd + 1;
        my $intronEnd = $start - 1;
        if($strandIndex eq $unStrn){
            push @{$blocks{$chrom}{$strandIndex}{$intronStart}{$intronEnd}}, ['intron', $strand, $name, $spliceIndex];
        } else {
            my $antisenseStrandIndex = $strandIndex eq '+' ? '-' : '+';
            push @{$blocks{$chrom}{$strandIndex}         {$intronStart}{$intronEnd}}, ['intron_sense', $strandIndex, $name, $spliceIndex];
            push @{$blocks{$chrom}{$antisenseStrandIndex}{$intronStart}{$intronEnd}}, ['intron_antisense', $antisenseStrandIndex, $name, $spliceIndex];
        }
    } 
}
#===============================================================================

#===============================================================================
#order splice elements (exons and introns) by position and trim off any portions that overlap BETWEEN genes
#-------------------------------------------------------------------------------
sub parseSegments { 
    print "parsing and trimming splice segments\n";
    open my $strandedH, ">", $ENV{MAP_STR_BED} or die "could not open $ENV{MAP_STR_BED} for writing: $!\n";
    open my $unstrandedH, ">", $ENV{MAP_UNSTR_BED} or die "could not open $ENV{MAP_UNSTR_BED} for writing: $!\n";
    foreach my $chrom(sort {$a cmp $b} keys %blocks){
        foreach my $strandIndex(sort {$a cmp $b} keys %{$blocks{$chrom}}){ 
            getSegments($chrom, $strandIndex, \my@segments);
            my ($prevStart, $prevEnd, @prevSegs) = (undef, undef);
            my $fileH = $strandIndex eq $unStrn ? $unstrandedH : $strandedH;
            my $bedStrand = $strandIndex eq $unStrn ? '+' : $strandIndex;
            $overlapBlock = ['ambiguous', $bedStrand, 'na', 0];
            foreach my $segment(@segments){
                my($start, $end) = @$segment;
                if(scalar @prevSegs){
                    if($start > $prevEnd){
                        commitSegments($fileH, $chrom, \@prevSegs);
                        @prevSegs = ($segment);
                        $prevStart = undef;
                        $prevEnd = undef;
                    } else { #overlap condition exists             
                        push @prevSegs, $segment;
                    }
                } else {
                    @prevSegs = ($segment);
                }
                $prevStart = ($prevStart and $prevStart < $start) ? $prevStart : $start;
                $prevEnd = ($prevEnd and $prevEnd > $end) ? $prevEnd : $end;
            }
            commitSegments($fileH, $chrom, \@prevSegs);
        }  
    }
    close $strandedH;
    close $unstrandedH;
}
sub getSegments { #sort by position the candidate segments on a chromosome strand
    my ($chrom, $strandIndex, $segments) = @_;
    foreach my $start(sort {$a <=> $b} keys %{$blocks{$chrom}{$strandIndex}}){
        foreach my $end(sort {$a <=> $b} keys %{$blocks{$chrom}{$strandIndex}{$start}}){
            foreach my $block(@{$blocks{$chrom}{$strandIndex}{$start}{$end}}){
                push @$segments, [$start, $end, $block];
            }
        }
    }
}
sub commitSegments {
    my ($fileH, $chrom, $segments) = @_;
    scalar(@$segments) > 1 and $segments = parseOverlap($segments);
    foreach my $segment(@$segments){ commitSegment($fileH, $chrom, $segment) }
}
sub parseOverlap { #when candidate segments overlap, divide into smaller segments while giving preference to sense strand elements
    my ($segments) = @_;
    splitOverlap($segments, \my@splitStarts, \my@splitEnds);
    my (@chunks, $chunkStart, $prevBlock, $prevSplitEnd);
    while(scalar(@splitStarts)){
        my $splitStart = shift @splitStarts;
        my $splitEnd = shift @splitEnds;
        getSplitBlocks($segments, $splitStart, $splitEnd, \my%senseTypeCounts, \my%elementTypeCounts, \my%blocks);
        my $senseType =   getPreferredType(\%senseTypeCounts, 'sense', 'antisense');
        my $elementType = getPreferredType($elementTypeCounts{$senseType}, 'exon', 'intron');
        my $n = $elementTypeCounts{$senseType}{$elementType};
        my $block = $n == 1 ? $blocks{$senseType}{$elementType} : $overlapBlock;
        if($prevBlock){
            unless($block eq $prevBlock){
                push @chunks, [$chunkStart, $prevSplitEnd, $prevBlock];
                $chunkStart = $splitStart;
            } 
        } else {
            $chunkStart = $splitStart;
        }
        $prevBlock = $block;
        $prevSplitEnd = $splitEnd;
    } 
    $prevBlock and push @chunks, [$chunkStart, $prevSplitEnd, $prevBlock];
    return \@chunks;
}
sub splitOverlap {
    my ($segments, $splitStarts, $splitEnds) = @_;
    my (%splitStarts, %splitEnds);
    foreach my $segment(@$segments){
        my ($segmentStart, $segmentEnd) = @$segment;    
        defined $segmentStart or next;  
        $segmentEnd - $segmentStart >= 1E6 and next; #disallow large elements (introns) from overlaps
        $splitStarts{$segmentStart}++;
        $splitStarts{$segmentEnd + 1}++;
        $splitEnds{$segmentEnd}++;
        $splitEnds{$segmentStart - 1}++;
    }
    @$splitStarts = sort {$a <=> $b} keys %splitStarts;
    @$splitEnds =   sort {$a <=> $b} keys %splitEnds;
    pop @$splitStarts;
    shift @$splitEnds;
}
sub getSplitBlocks {
    my ($segments, $splitStart, $splitEnd, $senseTypeCounts, $elementTypeCounts, $blocks) = @_;
    foreach my $segment(@$segments){
        my ($segmentStart, $segmentEnd, $block) = @$segment;
        defined $segmentStart or next;            
        if($splitStart <= $segmentEnd and $segmentStart <= $splitEnd){
            my ($featureType) = @$block;
            my $senseType =   $featureType =~ m/antisense/ ? 'antisense' : 'sense';
            my $elementType = $featureType =~ m/exon/      ? 'exon'      : 'intron';
            $$senseTypeCounts{$senseType}++;
            $$elementTypeCounts{$senseType}{$elementType}++;
            $$blocks{$senseType}{$elementType} = $block;
        }
    }
}  
sub getPreferredType {
    my($counts, $preferred, $nonPreferred) = @_;
    return $$counts{$preferred} ? $preferred : $nonPreferred;
}
sub commitSegment { #finalize the segment and print
    my ($fileH, $chrom, $segment) = @_;
    my ($start, $end, $block) = @$segment;
    defined $start or return;  
    $start <= $end or return;    
    my ($featureType, $strand, $name, $spliceIndex) = @$block;
    my $geneLength = $name eq 'na' ? 0 : $geneLengths{$chrom}{$name};
    $geneLength or $geneLength = 0;
    print $fileH join("\t", $chrom, $start - 1, $end, $name, 0, $strand,
                            $geneLength, $featureType, $spliceIndex), "\n";                     
}  
#===============================================================================

