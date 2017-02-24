#!/bin/bash

#q    require $GENOME $UCSC_FILE $NAME $CHROM_COLUMN $START_COLUMN $END_COLUMN $BED_FILE
#q    option $NAME_COLUMN[] $SCORE_COLUMN[] $STRAND_COLUMN[]

#$    -N  $NAME\_$GENOME
#$    -wd $GENOME_DIR
#$    -l  vf=2G
#$    -l  h_rt=2:00:00

#PBS  -N  $NAME\_$GENOME
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=2gb
#PBS  -l  walltime=2:00:00

if [ -f $UCSC_FILE ]; then

    echo "creating BED file of $GENOME $NAME"  
    
    UCSC_FILE="`echo $UCSC_FILE | sed -e 's|.gz$||' -e 's|.zip$||'`"    
    slurp $UCSC_FILE | 
    perl -e '
        $chc = "'$CHROM_COLUMN'";  
        $stc = "'$START_COLUMN'";  
        $enc = "'$END_COLUMN'";  
        $nac = "'$NAME_COLUMN'";  
        $scc = "'$SCORE_COLUMN'";  
        $src = "'$STRAND_COLUMN'"; 
        while(<>){
            chomp;
            @f = split("\t");
            $na = $nac ? $f[$nac-1] : "*";
            $sc = $scc ? $f[$scc-1] : 0;
            $sr = $src ? $f[$src-1] : "+";
            print join("\t", $f[$chc-1], $f[$stc-1], $f[$enc-1], $na, $sc, $sr), "\n";
        }
    ' | 
    slurp -o $BED_FILE
    checkPipe
    snipFile $BED_FILE

    echo "done"
    
else
    echo "does not exist: $UCSC_FILE"
fi
    
