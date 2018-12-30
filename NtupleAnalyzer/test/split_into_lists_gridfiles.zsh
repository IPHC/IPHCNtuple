#!/bin/env zsh
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

#version="ttH2017"
version="tHq2017_v1"

fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleProducer/"$version"/"

dataStr="SingleElectron|SingleMuon|DoubleEG|DoubleMuon|MuonEG"

liDATA=($(/usr/bin/rfdir ${fpath} | egrep -e "$dataStr" | awk '{print $9}'))
liMC=($(/usr/bin/rfdir ${fpath} | egrep -v "$dataStr" | awk '{print $9}'))

fpathDATAXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")
fpathMCXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")

#--- Files per job
nFilesDATA=10
nFilesMC=5
outDir="lists/"

rm -rf ${outDir}
mkdir ${outDir}

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
  
  split -a 5 -l ${nFilesDATA} -d /tmp/tempDATA.txt /tmp/${d1}_
  lsfi=($(ls /tmp/${d1}_*))
  jid=0
  
  for fil in $lsfi
  do
    if [[ $#d2 != 1 ]]; then
      mv ${fil} ${outDir}${d1}_${id2}_ID${jid}.txt
    else
      mv ${fil} ${outDir}${d1}_ID${jid}.txt
    fi
     
    jid=$[$jid+1]
  done
   
  rm -f /tmp/tempDATA.txt
done


# --- MC ---

rm -f /tmp/tempMC.txt
for line in $liMC
do
  #TTbar_Single files are very large --> split further
  if [[ $line == TTJets_SingleLeptFromT*_ext* ]]; then
    nFilesMC_tmp=1
  else 
    nFilesMC_tmp=$nFilesMC
  fi

  echo $line
  
  d1=$(echo $line)
  liMC2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}')
  echo $liMC2 | while read line2
  do
        f1=$(echo $line2)
        file=$(echo ${fpathMCXRD}${d1}/${f1})
        echo "${file}" >> /tmp/tempMC.txt
  done
  
  split -a 5 -l ${nFilesMC_tmp} -d /tmp/tempMC.txt /tmp/${d1}_
  lsfi=($(ls /tmp/${d1}_*))
  jid=0
  
  for fil in $lsfi
  do
    sampStrip=$(echo $id2 | sed "s%_RunIISpring15MiniAODv2_.*%%g")
    if [[ $#d2 != 1 ]]; then
      mv ${fil} ${outDir}${d1}_${sampStrip}_ID${jid}.txt
    else
    	echo ${fil}
	
	echo ${outDir}${d1}_{sampStrip}
      mv ${fil} ${outDir}${d1}_ID${jid}.txt
    fi
    
    jid=$[$jid+1]
  done
  
  rm -f /tmp/tempMC.txt
done

