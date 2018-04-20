#!/bin/bash

#eosdir="/store/user/chanon/TZQ/TestNtuplesV7/"
dir=/tmp/chanon

#while read line
for inputfile in `cat FileList_Syst.txt`
do

  #echo $line > tmp
  #proc=`awk '{print $1}' tmp`
  #inputfile=`awk '{print $2}' tmp`
  #echo $proc $inputfile

  root -l -q 'SystMerger.C("'${dir}/${inputfile}'")' #| awk 'NR==3' > tmp
  #nEntries=`cat tmp`
  #echo nEntries=$nEntries

  #nEntriesTot=$(( nEntriesTot + nEntries ));
  #echo nEntriesTot=$nEntriesTot

done 

#cd /tmp/chanon/
#ls *.root > l

#for file in `cat l`
#do
#  cmsStage $file ${eosdir}/$file
#  echo $file
#done

#cd -

#echo nEntriesTot=$nEntriesTot
