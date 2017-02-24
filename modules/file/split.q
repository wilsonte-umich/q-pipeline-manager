
invoke file/require.q
invoke file/check_prefix.q 
invoke math/integer.q  $VALUE $N_LINES
invoke math/non_zero.q $VALUE $N_LINES
invoke file/gzipped.q

qsub split.sh

<file name=split.sh>
#!/bin/bash

#q require $FILE $N_LINES $PREFIX $GZIPPED
#q option  $RM_FILE[] $COMPRESS[]

#$    -N  split_$N_LINES
#$    -l  vf=2G
#$    -l  h_rt=4:00:00

#PBS  -N  split_$N_LINES
#PBS  -l  mem=2gb
#PBS  -l  walltime=4:00:00

echo "splitting $FILE into $PREFIX, $N_LINES lines per file"
if [ "$GZIPPED" = "" ]; then
    GUNZIP="cat"
else 
    GUNZIP="gunzip -c"
fi

echo 
echo $FILE
cat $FILE | $GUNZIP | head
echo "..."
echo 

echo "removing any prior split files"
rm $PREFIX.[0-9][0-9][0-9] 2>/dev/null
rm $PREFIX.[0-9][0-9][0-9].gz 2>/dev/null

echo "splitting"
slurp $FILE | $GUNZIP | split -a 3 -d -l $N_LINES - $PREFIX.
checkPipe

if [ "$COMPRESS" != "" ]; then
    echo "compressing split files"
    gzip $PREFIX.[0-9][0-9][0-9]
    checkPipe
fi

if [ "$RM_FILE" != "" ]; then
    echo "deleting $FILE"
    rm $FILE
fi
echo "done"

</file>

