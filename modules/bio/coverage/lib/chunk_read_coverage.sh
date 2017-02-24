#!/bin/bash
#q    require $NAME $CHUNKS $CHUNK_SIZES

#$    -N  chunk_readCov_$NAME
#$    -l  vf=4G
#$    -l  h_rt=3:00:00
#$    -t  1-$N_CHUNKS

#PBS  -N  chunk_readCov_$NAME
#PBS  -l  mem=4gb
#PBS  -l  walltime=3:00:00
#PBS  -t  1-$N_CHUNKS

checkTaskID
getTaskObject CHUNK $CHUNKS
getTaskObject CHUNK_SIZE $CHUNK_SIZES

echo "extracting read coverage for $BAM_FILE $CHUNK ($CHUNK_SIZE bases)"

OUT_FILE=$BAM_FILE.readCov_$CHUNK

samtools mpileup -r $CHUNK $BAM_FILE | 
awk '{
  cs[$4]++;
  bc++;
}END{
  cs[0]='$CHUNK_SIZE'-bc;
  for(c in cs){
    print c,cs[c];
  }
}' > $OUT_FILE
checkPipe
echo "done"

