#!/bin/bash

#q    require $GENOME_DIR $GENOME $FILES $URL $TYPE

#$    -N  get_$GENOME\_$TYPE
#$    -wd $GENOME_DIR
#$    -l  vf=1G
#$    -l  h_rt=4:00:00

#PBS  -N  get_$GENOME\_$TYPE
#PBS  -d  $GENOME_DIR
#PBS  -l  mem=1gb
#PBS  -l  walltime=4:00:00

echo "getting $GENOME $TYPE files"
echo

cd $GENOME_DIR
IFS=" "
for FILE in $FILES
do
    WGET="$URL/$FILE"
    wget --timestamping --no-verbose $WGET
    if [ -f $FILE ]; then
        chmod ug+rw $FILE
        if [[ $FILE == *.tar.gz ]]; then
            tar -xzf $FILE
        elif [[ $FILE == *.zip ]]; then
            unzip -q $FILE
        elif [[ $FILE == *.gz ]]; then
            OUT_FILE="`echo $FILE | sed 's|.gz$||'`"
            gunzip -c $FILE > $OUT_FILE
            chmod ug+rw $OUT_FILE
        fi
    fi
    sleep 5
done

echo
echo "done"

