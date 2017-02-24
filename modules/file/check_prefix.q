
$TEST_FILE $PREFIX.__Q__CHECK__PREFIX__
$ERROR run [ "`touch $TEST_FILE 2>&1; rm $TEST_FILE 2>&1`" = "" ] || \
  echo '"$__Q__MASTER__FILE__: invalid prefix: $PREFIX"'
dieIf $ERROR

