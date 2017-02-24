
invoke file/require.q

preserve $GZIPPED run [ "`file $FILE | grep gzip`" = "" ] || echo "TRUE"

