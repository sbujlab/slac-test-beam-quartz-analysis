
if [ $# -gt 0 ] ; then
	run_list=$*
else
	run_list=$(more values_by_run.csv | awk -F',' '{print $1}')
fi

for run_num in ${run_list} ; do
	detector=$(grep ${run_num}, values_by_run.csv | grep -v ,${run_num}, | awk -F',' '{print $4}')
	if [[ ${detector} == "tandem" ]] ; then
echo		root -l "pmt_analyzer_tandem.c(${run_num})"
	elif [[ ${detector} == "sam" ]] ; then
		continue
	elif [[ ${detector:3:4} == "s" ]] ; then 
		root -l "pmt_analyzer_stack.c(${run_num})"
	else
echo		root -l "pmt_analyzer.c(${run_num})"
	fi
	sleep 0.5
done
