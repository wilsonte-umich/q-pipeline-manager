#!/bin/bash
#q    require $BAM_FILE

#$    -N  merge_readCov_$NAME
#$    -l  vf=1G
#$    -l  h_rt=1:00:00

#PBS  -N  merge_readCov_$NAME
#PBS  -l  mem=1gb
#PBS  -l  walltime=1:00:00

echo "merging read coverages for $BAM_FILE"

FILE_GLOB=$BAM_FILE.readCov_*

cat $FILE_GLOB |
awk '{
  cs[$1]+=$2;
}END{
  for(c in cs){
    print c,cs[c];
  }
}' | 
sort -k1,1n
checkPipe

rm $FILE_GLOB

echo "done"


