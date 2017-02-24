


$IN_STREAM



$ARRAY_VALUES run
slurp $IN_STREAM | 
cut -f $ARRAY_COLUMN |
sort |
groupBy -g 1 -c 1 -o count |
cut -f 1 

$N_JOBS run echo $ARRAY_VALUES | wc -w

$ARRAY_COLUMNS
$ARRAY_VALUES





invoke $TARGET

-t $N_JOBS


