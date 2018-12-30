#!/bin/env zsh
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

#fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/FlatTree/tHq2017_v0/"
fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/FlatTree/tHq2017_v1/"

dataStr="SingleElectron|SingleMuon|DoubleEG|DoubleMuon|MuonEG"

liDATA=($(/usr/bin/rfdir ${fpath} | egrep -e "$dataStr" | awk '{print $9}'))
liMC=($(/usr/bin/rfdir ${fpath} | egrep -v "$dataStr" | awk '{print $9}'))

fpathDATAXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")
fpathMCXRD=$(echo ${fpath} | sed "s%/dpm%root://sbgse1.in2p3.fr//dpm%g")

nFilesDATA=40
nFilesMC=10
outDir="lists/"

rm -rf ${outDir}
mkdir ${outDir}

# DATA

rm -f /tmp/tempDATA.txt
for line in $liDATA
do
  echo $line
  d1=$(echo $line)
  liDATA2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}')
  d2=($(echo $liDATA2))
  for id2 in $d2
  do
    liDATA3=$(/usr/bin/rfdir ${fpath}${d1}/${id2})
    d3=$(echo $liDATA3 | awk '{print $9}')
    liDATA4=($(/usr/bin/rfdir ${fpath}${d1}/${id2}/${d3} | awk '{print $9}'))
    for d4 in $liDATA4
    do
      liDATA5=($(/usr/bin/rfdir ${fpath}${d1}/${id2}/${d3}/${d4} | awk '{print $9}'))
      for line2 in $liDATA5
      do
        f1=$(echo $line2)
        file=$(echo ${fpathDATAXRD}${d1}/${id2}/${d3}/${d4}/${f1})
        echo "${file}" >> /tmp/tempDATA.txt
      done
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
done

# MC

rm -f /tmp/tempMC.txt
for line in $liMC
do
  echo $line
  d1=$(echo $line)
  liMC2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}')
  d2=($(echo $liMC2))
  for id2 in $d2
  do
    liMC3=$(/usr/bin/rfdir ${fpath}${d1}/${id2})
    d3=$(echo $liMC3 | awk '{print $9}')
    liMC4=($(/usr/bin/rfdir ${fpath}${d1}/${id2}/${d3} | awk '{print $9}'))
    for d4 in $liMC4
    do
      liMC5=($(/usr/bin/rfdir ${fpath}${d1}/${id2}/${d3}/${d4} | awk '{print $9}'))
      for line2 in $liMC5
      do
        f1=$(echo $line2)
        file=$(echo ${fpathMCXRD}${d1}/${id2}/${d3}/${d4}/${f1})
        echo "${file}" >> /tmp/tempMC.txt
      done
    done
    split -a 5 -l ${nFilesMC} -d /tmp/tempMC.txt /tmp/${d1}_
    lsfi=($(ls /tmp/${d1}_*))
    jid=0
    for fil in $lsfi
    do
      sampStrip=$(echo $id2 | sed "s%_RunIISpring15MiniAODv2_.*%%g")
      if [[ $#d2 != 1 ]]; then
        mv ${fil} ${outDir}${d1}_${sampStrip}_ID${jid}.txt
      else
          mv ${fil} ${outDir}${d1}_ID${jid}.txt
      fi
      jid=$[$jid+1]
    done
    
    rm -f /tmp/tempMC.txt
  done
done


### REMOVE MANUALLY SOME SAMPLES HERE (keep FlatTree files, but don't run NTProd on them bc not needed)
rm ./lists/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8*
rm ./lists/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8*
rm ./lists/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8*
rm ./lists/TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8*
rm ./lists/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8*
