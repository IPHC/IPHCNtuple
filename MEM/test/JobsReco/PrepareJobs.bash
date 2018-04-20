#!/bin/bash

#Automatically creates jobs, given the processes and files listed in FileList.txt
#Environment variables LHAPDF and MEMEXECDIR (containing the executable) must be set

#Suffix of the subdirectory: jobs will be created in Jobs_${opt}, that should contain the config.cfg file 
opt=thq_test

#nEv events are run per job. Recommended nEv=6. If running also TTWJJ hyp, use nEv=1 (very slow). Depends on process
nEv=20

#LSF queue
queue=1nd

mkdir Jobs_${opt}
cp config.cfg Jobs_${opt}/config.cfg

while read line
do

  echo $line > tmp
  proc=`awk '{print $1}' tmp`
  inputfile=`awk '{print $2}' tmp`
  echo $proc $inputfile

  root -l -q 'ReadEntries.C("'${inputfile}'")' | awk 'NR==3' > tmp
  nEntries=`cat tmp`

  nJobs=$(($nEntries / $nEv ))
  echo nEntries=$nEntries nJobs=$nJobs

  #nJobsMax=$((20000/$nEv))
  #echo nJobsMax=$nJobsMax

  cp RunBatchMEM_template Jobs_${opt}/RunBatchMEM_${proc}_template
  sed -i s/OPTION/${opt}/g Jobs_${opt}/RunBatchMEM_${proc}_template
  sed -i s/PROC/${proc}/g Jobs_${opt}/RunBatchMEM_${proc}_template
  sed -i s,INPUTFILE,${inputfile},g Jobs_${opt}/RunBatchMEM_${proc}_template
  sed -i s,DIR_JOBS,`pwd`,g Jobs_${opt}/RunBatchMEM_${proc}_template
  sed -i s,DIR_LHAPDF,${LHAPDF},g Jobs_${opt}/RunBatchMEM_${proc}_template
  sed -i s,DIR_MEMEXEC,${MEMEXECDIR},g Jobs_${opt}/RunBatchMEM_${proc}_template #MEMEXECDIR must be defined!

  for i in `seq 0 $nJobs` #$nJobs 
  do
#    if [[ $i -ge $nJobsMax ]]; then break; fi

    i1=$(($i*$nEv))
    i2=$(($i*$nEv+$nEv))
    cp Jobs_${opt}/RunBatchMEM_${proc}_template Jobs_${opt}/RunBatchMEM_${proc}_$i 
    sed -i s/NMIN/$i1/g Jobs_${opt}/RunBatchMEM_${proc}_$i
    sed -i s/NMAX/$i2/g Jobs_${opt}/RunBatchMEM_${proc}_$i
    sed -i s/NUM/$i/g Jobs_${opt}/RunBatchMEM_${proc}_$i
    echo bsub -q ${queue} -N RunBatchMEM_${proc}_$i
  done
done < FileList.txt

chmod +x Jobs_${opt}/RunBatchMEM_*
