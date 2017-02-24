
invoke bio/genome/file_schema.q  
$STAT_VALUE run get_stat $GENOME_DIR $GENOME $STAT_NAME

preserve $$STAT_NAME $STAT_VALUE

