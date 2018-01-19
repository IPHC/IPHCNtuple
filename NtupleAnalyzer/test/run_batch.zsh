#!/bin/env zsh

# don't forget /opt/sbg/scratch1/cms

echo "Don't forget to update the lumi and the maximum number of events to run on in this script if needed !"

isdata=0
doSystCombine=0
lumi=35900

cp /tmp/x509up_u8066 /home-pbs/ntonon/proxy

jName=${1}
if [[ ${jName} == "" ]]; then
  echo "Please specify the run name"
  exit 1
fi

#que="cms_local"
que="cms"

export HOME=$(pwd)

dout="/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleAnalyzer/test"
dout_f="/opt/sbg/scratch1/cms/ntonon/Analyzer_ntuples_prod_walrus_patch2/"

runName="toy_${jName}"
logName="log_${jName}"

rm -rf ${dout_f}/${runName}
mkdir ${dout_f}/${runName}
rm -rf ${logName}
mkdir ${logName}

nmax=-1

fxsec="table.txt"

#fdir=$(ls -d lists_MC_Walrus-patch1)
fdir=$(ls -d lists_tHq)

echo $fdir

echo $fdir | while read line
do
fpath="${HOME}/${line}/"
flist=$(ls ${fpath})
dir=${line}

echo $flist | while read line
do
  jidx=0
  sample=$(echo $line | sed 's%.txt%%g')
  dataset=$(echo $sample | sed 's%_ID..*%%g')
  if [[ ! -d ${runName}/${dataset} ]]; then
    mkdir ${dout_f}/${runName}/${dataset}
  fi
  linexsec=$(grep $dataset $fxsec)
  nowe=$(echo $linexsec | awk '{print $3}')
  xsec=$(echo $linexsec | awk '{printf $2}')
  if [[ $nowe == "" ]]; then
    nowe=1
  fi
  if [[ $xsec == "" ]]; then
    xsec=1
  fi
  fl=$(echo $sample | cut -c1-1)
  fl2=$(echo $sample | cut -c1-11)
  datamc=""
  if [[ $fl == "Run2016" ]]; then
    isdata=1
    doSystCombine=0
    datamc="DATA"
    nmax=${nmax}
  else
    isdata=0
    doSystCombine=0
    datamc="MC"
    nmax=${nmax}
  fi
  
  
  fout=`echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g'`
  lout=`echo ${line}_${jidx} | sed 's%.txt%%g'`
  #fout=$(echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g')
  #lout=$(echo ${line}_${jidx} | sed 's%.txt%%g')

  echo "${dataset}: $nowe $xsec $lumi"
  #echo "${fpath}${line}"
 
  qsub -N ${dir} -q ${que} -o ${logName}/${sample}.log -j oe single_batch_job.sh \
-v dout=${dout},line2=${fpath}${line},dout_f=${dout_f},fout=${fout},nowe=${nowe},xsec=${xsec},lumi=${lumi},isdata=${isdata},doSystCombine=${doSystCombine},dataset=${dataset},nmax=${nmax}
done

echo "going to sleep 2700 s (45 mn)"
sleep 2700

done 
