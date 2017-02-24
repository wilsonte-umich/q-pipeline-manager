
$FILE $READ_GLOB
invoke file/require.q
invoke file/gzipped.q
invoke file/create.q $DIR $OUTPUT_DIR

qsub fastqc.sh

<file name=fastqc.sh>
#!/bin/bash

#q    require $NAME $READ_GLOB $OUTPUT_DIR $GZIPPED
#q    option  $TMP_DIR[/tmp] 

#$    -N  fastqc.$NAME
#$    -wd $OUTPUT_DIR
#$    -l  vf=2G
#$    -l  h_rt=2:00:00

#PBS  -N  fastqc.$NAME
#PBS  -d  $OUTPUT_DIR
#PBS  -l  mem=2gb
#PBS  -l  walltime=2:00:00

echo "running FastQC on $NAME"

echo "creating tmp file"
mkdir -p $TMP_DIR
FILE_NAME="$OUTPUT_DIR/$NAME.fastqc"
FILE_NAME="`echo $FILE_NAME | sed 's|/|_|g'`"
TMP_FILE="$TMP_DIR/$FILE_NAME"
if [ "$GZIPPED" != "" ]; then
    slurp -o $TMP_FILE gunzip -c $READ_GLOB
else
    slurp $READ_GLOB | slurp -o $TMP_FILE
fi
checkPipe

echo "running FastQC"
echo
fastqc --outdir $OUTPUT_DIR $TMP_FILE
echo

echo "cleaning up tmp file"
rm $TMP_FILE

echo "renaming output"
mv -f $OUTPUT_DIR/$FILE_NAME"_fastqc"     $OUTPUT_DIR/$NAME"_fastqc"
mv -f $OUTPUT_DIR/$FILE_NAME"_fastqc.zip" $OUTPUT_DIR/$NAME"_fastqc.zip"

echo "done"

</file>

