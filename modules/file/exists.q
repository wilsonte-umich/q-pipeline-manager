
$N_FILES run ls -1d $FILE 2>/dev/null | wc -l | sed 's/^0$//'
$EXISTS  $N_FILES ? TRUE : FALSE
$N_FILES $N_FILES ? $N_FILES : 0

preserve $EXISTS $EXISTS
preserve $N_FILES $N_FILES

