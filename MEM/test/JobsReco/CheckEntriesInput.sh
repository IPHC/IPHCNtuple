#!/bin/bash

nEntriesTot=0

while read line
do

  echo $line > tmp
  proc=`awk '{print $1}' tmp`
  inputfile=`awk '{print $2}' tmp`
  echo $proc $inputfile

  root -l -q 'ReadEntries.C("'${inputfile}'")' | awk 'NR==3' > tmp
  nEntries=`cat tmp`
  echo nEntries=$nEntries

  nEntriesTot=$(( nEntriesTot + nEntries ));
  #echo nEntriesTot=$nEntriesTot

done < FileList.txt

echo nEntriesTot=$nEntriesTot
