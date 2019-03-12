

for run_num in $(more values_by_run.csv | awk -F',' '{print $1}') ; do
	signal=$(grep ${run_num}, values_by_run.csv | awk -F',' '{print $7}')
	if [ ${signal} -gt 300 ] ; then
		root -l "pmt_analyzer.c(${run_num})"
		sleep 1
	fi
done
