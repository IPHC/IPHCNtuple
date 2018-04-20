#!/bin/bash

#Suffix of the subdirectory: jobs will be created in Jobs_${opt}, that should contain the config.cfg file 
opt=thq_test

#nEv events are run per job. Recommended nEv=6. If running also TTWJJ hyp, use nEv=1 (very slow).
nEv=20

#LSF queue
queue=1nd

minimumsize=5000

while read line
do

  echo $line > tmp
  proc=`awk '{print $1}' tmp`
  inputfile=`awk '{print $2}' tmp`
#  echo $proc $inputfile

  root -l -q 'ReadEntries.C("'${inputfile}'")' | awk 'NR==3' > tmp
  nEntries=`cat tmp`

  #nJobs=2000
  nJobs=$(($nEntries / $nEv ))
#  echo nEntries=$nEntries nJobs=$nJobs

  for i in `seq 0 $nJobs` #$nJobs 
  do
    #if [ ! -s Jobs_${opt}/output_${proc}_${opt}_${i}.root ] ; then rm -f Jobs_${opt}/output_${proc}_${opt}_${i}.root ; fi
    actualsize=$(wc -c <"Jobs_${opt}/output_${proc}_${opt}_${i}.root") 
    #echo Jobs_${opt}/output_${proc}_${opt}_${i}.root $actualsize
    if [  $actualsize -le $minimumsize ] ; then rm Jobs_${opt}/output_${proc}_${opt}_${i}.root ; fi
    if [ ! -e Jobs_${opt}/output_${proc}_${opt}_${i}.root ] ; then echo bsub -q ${queue} RunBatchMEM_${proc}_${i} ; fi
  done

done < FileList.txt



