#NB : WILL COMPILE CODE WHEN JOB STARTS TO RUN ! SO IF E.G. YOU CHANGE THE
#SAMPLES TO PROCESS IN BETWEEN, YOU WILL GET DIFFERENT RESULT THAN EXPECTED !!

#NB2 : use cms_local_mdm for <4h jobs, cms_local for <72h jobs

#NB3 : must give a name (for logfile) as argument

if [[ $1 == "" ]]; then
	echo "Missing arg ! Exit"
	exit 1
fi

cd ..
make
cd test
qsub -q cms_local -o log_PU/log_$1.log launch_job_PU.sh
