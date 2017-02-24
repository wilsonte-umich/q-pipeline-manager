
invoke file/exists.q
preserve $N_FILES $N_FILES
$NOT_REQUIRED $NO_REQUIRE ? \
  run IFS=" "; for f in "$NO_REQUIRE"; do [ "\$f" = "$FILE" ] && echo 1; done : \
  $NULL
exitIf $NOT_REQUIRED
$ERROR $EXISTS ? FALSE : "$__Q__MASTER__FILE__: no matching files found: $FILE" 
dieIf $ERROR

