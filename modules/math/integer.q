
$ERROR run [ "`echo $VALUE | tr -d [:digit:]`" = "" ] || \
  echo '"$__Q__MASTER__FILE__: $VALUE is not an integer"'
dieIf $ERROR

