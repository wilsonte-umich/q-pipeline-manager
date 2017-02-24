
# set the variables and paths
invoke lib/common.q
$MIN_MAPQ   $MIN_MAPQ ? $MIN_MAPQ : 10;

# find the coverage breaks for each sample
qsub lib/make_map.sh
