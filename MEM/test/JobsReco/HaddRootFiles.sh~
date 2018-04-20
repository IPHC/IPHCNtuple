#!/bin/bash

#Suffix of the subdirectory: jobs will be created in Jobs_${opt}, that should contain the config.cfg file 
opt=thq_test

while read line
do

  echo $line > tmp
  proc=`awk '{print $1}' tmp`
#  inputfile=`awk '{print $2}' tmp`
#  echo $proc $inputfile

#  root -l -q 'ReadEntries.C("'${inputfile}'")' | awk 'NR==3' > tmp
#  nEntries=`cat tmp`

#  nJobs=$(($nEntries / $nEv ))
#  echo nEntries=$nEntries nJobs=$nJobs

  hadd Jobs_${opt}/output_${proc}_${opt}_all.root Jobs_${opt}/output_${proc}_${opt}_*.root

#  for i in `seq 0 $nJobs` #$nJobs 
#  do
#    if [ ! -e Jobs_${opt}/output_${proc}_${opt}_${i}.root ] ; then echo bsub -q ${queue} RunBatchMEM_${proc}_${i} ; fi
#  done

done < FileList.txt
