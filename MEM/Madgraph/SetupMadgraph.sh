#!/bin/bash

ls -d */ | awk -F "/" '{ print $1 }' > ProcessList.txt

for dir in `cat ProcessList.txt`
do

  cd ${dir}/src
  make clean
  make

  cd ../SubProcesses
  ls -d */ | awk -F "/" '{ print $1 }' > listdir
  for dir in `cat listdir`
  do
    cd ${dir}
    make clean
    make
    cd ..
  done
  cd ../..

done

