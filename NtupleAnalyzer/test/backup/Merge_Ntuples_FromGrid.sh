### MERGE NTUPLE FILES (FROM GRID TO LOCAL) ###

#!/bin/env zsh
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

#Path of directories containing Ntuple files to merge
fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/" 

dataStr="SingleElectron|SingleMuon|DoubleEG|DoubleMuon|MuonEG"

#[awk '{print $9}']  --> list of all the filenames ((filename is 9th argument of ls -l))
#Class dirs. as MC or Data depending if they contrain dataStr or not
liDATA=($(/usr/bin/rfdir ${fpath} | egrep -e "$dataStr" | awk '{print $9}'))
liMC=$(/usr/bin/rfdir ${fpath} | egrep -v "$dataStr" | awk '{print $9}')

echo $liMC

#Add sbgse prefix to files
fpathDATAXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")
fpathMCXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")

#Create local dir. to store merged ntuples
outDir=$PWD/merged_ntuples
mkdir $outDir

# --- MC ---
rm -f /tmp/tempMC.txt

#List all MC dirs.
for line in $liMC 
do
  echo $line
  d1=$(echo $line)
  
  #List all filenames in dir.
  liMC2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}') 
  #echo $liMC2
  
  for line2 in $liMC2
  do
        f1=$(echo $line2)
	
	#Get complete filepath
        file=$(echo ${fpathMCXRD}${d1}/${f1}) 
	#echo $file
	 
	#Get full list of files into tmp file
        echo "${file}" >> /tmp/tempMC.txt 
  done

  #Write hadd command and full list of files into tmp script ; execute it and delete it
  echo "hadd -f $outDir/merged_$line.root " > /tmp/tmp_ScriptMerging_$line.sh
  cat /tmp/tempMC.txt >> /tmp/tmp_ScriptMerging_$line.sh
  tr '\n' ' ' < /tmp/tmp_ScriptMerging_$line.sh > /tmp/ScriptMerging_$line.sh
  rm -f /tmp/tmp_ScriptMerging_$line.sh
  chmod 755 /tmp/ScriptMerging_$line.sh
  #/tmp/ScriptMerging_$line.sh
  #rm -f /tmp/ScriptMerging_$line.sh
  
  rm -f /tmp/tempMC.txt
done

# --- DATA --- (idem)
rm -f /tmp/tempDATA.txt
for line in $liDATA
do
  echo $line
  d1=$(echo $line)
  liDATA2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}')

  for line2 in $liDATA2
  do
        f1=$(echo $line2)
        file=$(echo ${fpathDATAXRD}${d1}/${f1})
        echo "${file}" >> /tmp/tempDATA.txt
  done
  
  echo "hadd -f  $outDir/merged_$line.root " > /tmp/tmp_ScriptMerging_$line.sh
  cat /tmp/tempDATA.txt >> /tmp/tmp_ScriptMerging_$line.sh
  tr '\n' ' ' < /tmp/tmp_ScriptMerging_$line.sh > /tmp/ScriptMerging_$line.sh
  rm -f /tmp/tmp_ScriptMerging_$line.sh
  chmod 755 /tmp/ScriptMerging_$line.sh
  #/tmp/ScriptMerging_$line.sh
  #rm -f /tmp/ScriptMerging_$line.sh

   
  rm -f /tmp/tempDATA.txt
done


