#!/bin/env zsh
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066
fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/"

dataStr="SingleElectron|SingleMuon|DoubleEG|DoubleMuon|MuonEG"

liDATA=($(/usr/bin/rfdir ${fpath} | egrep -e "$dataStr" | awk '{print $9}'))
liMC=($(/usr/bin/rfdir ${fpath} | egrep -v "$dataStr" | awk '{print $9}'))

fpathDATAXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")
fpathMCXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")

outDir=$PWD/merged_ntuples
mkdir $outDir

# --- MC ---

rm -f /tmp/tempMC.txt
for line in $liMC
do
  echo $line
  d1=$(echo $line)
  liMC2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}') #lists all FILENAMES in dir
  echo $liMC2 | while read line2
  do
        f1=$(echo $line2)
        file=$(echo ${fpathMCXRD}${d1}/${f1})
        echo "${file}" >> /tmp/tempMC.txt
  done

  echo "hadd $outDir/merged_$line.root " > /tmp/tmp_ScriptMerging_$line.sh
  cat /tmp/tempMC.txt >> /tmp/tmp_ScriptMerging_$line.sh
  tr '\n' ' ' < /tmp/tmp_ScriptMerging_$line.sh > /tmp/ScriptMerging_$line.sh
  rm -f /tmp/tmp_ScriptMerging_$line.sh
  chmod 755 /tmp/ScriptMerging_$line.sh
  /tmp/ScriptMerging_$line.sh
  rm -f /tmp/ScriptMerging_$line.sh
  
  rm -f /tmp/tempMC.txt
done

# --- DATA ---

rm -f /tmp/tempDATA.txt
for line in $liDATA
do
  echo $line
  d1=$(echo $line)
  liDATA2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}')
  echo $liDATA2 | while read line2
  do
        f1=$(echo $line2)
        file=$(echo ${fpathDATAXRD}${d1}/${f1})
        echo "${file}" >> /tmp/tempDATA.txt
  done
  
  echo "hadd $outDir/merged_$line.root " > /tmp/tmp_ScriptMerging_$line.sh
  cat /tmp/tempDATA.txt >> /tmp/tmp_ScriptMerging_$line.sh
  tr '\n' ' ' < /tmp/tmp_ScriptMerging_$line.sh > /tmp/ScriptMerging_$line.sh
  rm -f /tmp/tmp_ScriptMerging_$line.sh
  chmod 755 /tmp/ScriptMerging_$line.sh
  /tmp/ScriptMerging_$line.sh
  rm -f /tmp/ScriptMerging_$line.sh

   
  rm -f /tmp/tempDATA.txt
done

hadd $outDir/SingleElectron.root $outDir/*SingleElectron*
hadd $outDir/SingleMuon.root $outDir/*SingleMuon*
hadd $outDir/MuonEG.root $outDir/*MuonEG*
hadd $outDir/DoubleEG.root $outDir/*DoubleEG*
hadd $outDir/DoubleMuon.root $outDir/*DoubleMuon*

rm $outDir/*SingleElectron*2016*.root
rm $outDir/*SingleMuon*2016*.root
rm $outDir/*DoubleEG*2016*.root
rm $outDir/*DoubleMuon*2016*.root
rm $outDir/*MuonEG*2016*.root
