#!/bin/env zsh

#version="ttH2017" #output subdir
version="tHq2017" #output subdir

echo ""
#Don't forget to update the lumi and the xsec table

isdata=0
doSystCombine=0
lumi=41500

cp /tmp/x509up_u8066 /home-pbs/ntonon/proxy

jName=${1}
if [[ ${jName} == "" ]]; then
  echo "Please specify the run name"
  exit 1
fi

#que="cms_local_mdm"
#que="cms"
que="cms_local_short" #reseved slots, fast, <1h jobs only


export HOME=$(pwd)

dout="/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test" #Current dir.
dout_f="/opt/sbg/scratch1/cms/ntonon/Analyzer_ntuples_tHq" #tmp output dir

runName="toy_${jName}"
logName="log_${jName}"

rm -rf ${dout_f}/${runName}
mkdir ${dout_f} #NEW
mkdir ${dout_f}/${runName}
rm -rf ${logName}
mkdir ${logName}

nmax=-1

fxsec="table_MC_tHqAnalysis.txt"

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
  dataset=$(echo $sample | sed 's%__ID..*%%g') #CHANGED
  if [[ ! -d ${runName}/${dataset} ]]; then
    mkdir ${dout_f}/${runName}/${dataset} 
  fi
  linexsec=$(grep -m 1 $dataset $fxsec) #grep -m 1 <-> only first occurence
  nowe=$(echo $linexsec | awk '{print $3}')
  xsec=$(echo $linexsec | awk '{printf $2}')
  
  if [[ $nowe == "" ]]; then
    nowe=1
  fi
  if [[ $xsec == "" ]]; then
    xsec=1
  fi
  datamc=""
  if [[ $sample == *Run2017* ]]; then
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
  
    
  #isdata=1
  
  fout=`echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g'`
  lout=`echo ${line}_${jidx} | sed 's%.txt%%g'`
  #fout=$(echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g')
  #lout=$(echo ${line}_${jidx} | sed 's%.txt%%g')

  echo "* ${dataset} :"
  echo "- nowe = $nowe"
  echo "- xsec = $xsec" 
  echo "- lumi = $lumi"
  echo "- isdata = " ${isdata}
  #echo "${fpath}${line}"
 
  qsub -N ${dir} -q ${que} -o ${logName}/${sample}.log -j oe single_batch_job.sh \
-v dout=${dout},line2=${fpath}${line},dout_f=${dout_f},fout=${fout},nowe=${nowe},xsec=${xsec},lumi=${lumi},isdata=${isdata},doSystCombine=${doSystCombine},dataset=${dataset},nmax=${nmax},version=${version}
done

echo "going to sleep 2700 s (45 mn)"
sleep 2700

done 
