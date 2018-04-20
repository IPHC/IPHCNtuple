#!/bin/bash

#PROC_SA_CPP_sm_DECAY_gqlnullgq
model=sm_no_b_mass
suffix="ppthq"
Name="THQ"

ProcDir=PROC_SA_CPP_${model}_DECAY_${suffix}
echo Generating code for ${ProcDir} with MakeFile name ${Name}
echo
echo \#${Name}
echo ${Name}_MADLIBDIR=../Madgraph/${ProcDir}/lib
echo ${Name}_MADINCDIR=../Madgraph/${ProcDir}/src
echo ${Name}_MADMODELLIB=model_${model}_${suffix}

cd ${ProcDir}/SubProcesses
ls -d */ | awk -F "/" '{ print $1 }' > l

NUM=0
for dir in `cat l`
do
  #echo $NUM 
  echo ${Name}_MADPROCESS_${NUM}=CPPProcess_${dir}
  echo ${Name}_MADPROCDIR_${NUM}=../Madgraph/${ProcDir}/SubProcesses/${dir}
  ((NUM+=1))
done

((NUM-=1))

echo
echo Add to CXXFLAGS:
for i in `seq 0 ${NUM}`
do
  echo -n ' '-I\${${Name}_MADPROCDIR_${i}} 
done
echo -n ' '-I\${${Name}_MADINCDIR}
echo

echo
echo Add to LIBS:
echo ' '-L\${${Name}_MADLIBDIR} -l\${${Name}_MADMODELLIB}

echo
echo Add to objects:
for i in `seq 0 ${NUM}`
do
echo -n ' '\${${Name}_MADPROCDIR_${i}}/\${${Name}_MADPROCESS_${i}}.o
done
echo

echo
echo Add to TestMEMcomb.o:
echo ' '\${${Name}_MADLIBDIR}/lib\${${Name}_MADMODELLIB}.a

echo
echo Add at the end:
echo
for i in `seq 0 ${NUM}`
do
echo \${${Name}_MADPROCDIR_${i}}/\${${Name}_MADPROCESS_${i}}.o: \${${Name}_MADPROCDIR_${i}}/\${${Name}_MADPROCESS_${i}}.cc \${${Name}_MADLIBDIR}/lib\${${Name}_MADMODELLIB}.a
echo -e '\t'\${CXX} -l\${LIBS} \${CXXFLAGS} -c \${${Name}_MADPROCDIR_${i}}/\${${Name}_MADPROCESS_${i}}.cc
echo
done

echo Add to MEPhaseSpace.h:
for dir in `cat l`
do
  echo \#include \"../Madgraph/${ProcDir}/SubProcesses/${dir}/CPPProcess_${dir}.h\"
done

echo
echo Add to the MEPhaseSpace class:
for dir in `cat l`
do
  echo CPPProcess_${dir}* process_${dir}\;
done

echo
echo Add to the MEPhaseSpace constructor:
for dir in `cat l`
do
  echo MGcard = MadgraphDir + \"/${ProcDir}/Cards/param_card.dat\"\;
  echo process_${dir} = new CPPProcess_${dir}\(\)\;
  echo process_${dir}-\>initProc\(MGcard.c_str\(\)\)\;
#  echo process_${dir}-\>initProc\(\"/afs/cern.ch/work/c/chanon/MEM/Madgraph/${ProcDir}/Cards/param_card.dat\"\)\;
  echo "if (verbosity>=1) cout" \<\< \"${Name} Process nexternal=\" \<\< process_${dir}-\>nexternal \<\< endl\;
  echo
done

echo Manually change ComputeSubMatrixElement and the phase space computation




 
