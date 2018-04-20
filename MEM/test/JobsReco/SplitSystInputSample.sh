#!/bin/bash

#eosdir="/store/user/chanon/TZQ/TestNtuplesV7_Syst"
dir=/tmp/chanon

#for inputfile in `cat FileList.txt`
while read line
do

  echo $line > tmp
  proc=`awk '{print $1}' tmp`
  inputfile=`awk '{print $2}' tmp`
  echo $proc $inputfile

  echo ${inputfile}
#  echo $line > tmp
#  proc=`awk '{print $1}' tmp`
#  inputfile=`awk '{print $2}' tmp`
#  echo $proc $inputfile

  root -l -q 'SplitSyst.C("'${inputfile}'")' #| awk 'NR==3' > tmp
  #nEntries=`cat tmp`
  #echo nEntries=$nEntries

  #nEntriesTot=$(( nEntriesTot + nEntries ));
  #echo nEntriesTot=$nEntriesTot

done < FileList.txt

#cd /tmp/chanon/
#ls *.root > l

#for file in `cat l`
#do
#  cmsStage $file ${eosdir}/$file
#  echo $file
#done

#cd -

#echo nEntriesTot=$nEntriesTot
