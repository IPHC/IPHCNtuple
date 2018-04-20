#!/bin/bash

opt=test_tZqAllSamplesMoriondNew2_FakesNew2_minus

eosdir="/store/user/chanon/TZQ/TestNtuplesV7_Syst_MEMoutput"

while read line
do

  echo $line > tmp
  proc=`awk '{print $1}' tmp`
  inputfile=`awk '{print $2}' tmp`
#  echo $proc $inputfile

  ls Jobs_${opt}/output_${proc}_${opt}_all.root
  cmsStage Jobs_${opt}/output_${proc}_${opt}_all.root ${eosdir}/FCNCNTuple_${proc}_.root
  echo ${eosdir}/FCNCNTuple_${proc}.root

done < FileList.txt
