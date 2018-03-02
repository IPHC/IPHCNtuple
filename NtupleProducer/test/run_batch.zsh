#!/bin/env zsh

#--- UPDATE proxy
cp /tmp/x509up_u8066 /home-pbs/ntonon/proxy/

jName=${1}
if [[ ${jName} == "" ]]; then
  echo "Please specify the run name"
  exit 1
fi

que="cms"
#que="cms_local_mdm" #local queue

export HOME=$(pwd)

dout="/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleProducer/test"

dout_f="/opt/sbg/scratch1/cms/ntonon/ntuples_prod_tHq"


echo "CMSSW_RELEASE_BASE" $CMSSW_RELEASE_BASE

runName="toy${jName}"
logName="log${jName}"

rm -rf ${logName}
mkdir ${logName}
rm -rf ${dout_f}/${runName}
mkdir ${dout_f}
mkdir ${dout_f}/${runName}

nmax=-1

#fxsec="table.txt"

fdir=$(ls -d lists_tHq*)

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
    #mkdir ${runName} #removed
    #mkdir ${runName}/${dataset}  
    mkdir ${dout_f}/${runName}/${dataset}
  fi
  #linexsec=$(grep $dataset $fxsec)
  #noe=$(echo $linexsec | awk '{print $3}')
  #xsec=$(echo $linexsec | awk '{print $2}')
  #if [[ $noe == "" ]]; then
    #noe=1
  #fi
  #if [[ $xsec == "" ]]; then
    #xsec=1
  #fi
  datamc=""
  if [[ $sample == *Run2016* ]]; then
    isdata=1
    datamc="DATA"
    nmax=${nmax}
  else
    isdata=0
    datamc="MC"
    nmax=${nmax}
  fi
  
  #isdata=1
   
  fout=$(echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g')
  lout=$(echo ${line}_${jidx} | sed 's%.txt%%g')

# echo "${dataset}: $noe $xsec"
  echo "${fpath}${line}"
  echo "isdata = " ${isdata}
 
 # qsub -N ${dir} -q ${que} -o ${logName}/${sample}.log -j oe single_batch_job.sh \
#-v dout=${dout},line2=${fpath}${line},fout=${fout},noe=${noe},xsec=${xsec},isdata=${isdata},sample=${sample},nmax=${nmax},dout_f=${dout_f}


  qsub -N ${dir} -q ${que} -o ${logName}/${sample}.log -j oe single_batch_job.sh \
-v dout=${dout},line2=${fpath}${line},fout=${fout},isdata=${isdata},sample=${sample},nmax=${nmax},dout_f=${dout_f},dataset=${dataset}
done

#echo ${dataset}

echo "going to sleep 2700 s (45 mn)"
sleep 2700

done 
