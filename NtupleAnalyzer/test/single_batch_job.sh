#!/bin/sh

export X509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/
export SCRAM_ARCH=slc6_amd64_gcc630
eval `scramv1 runtime -sh`
cd -

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${dout}/../:${dout}/../../NtupleProducer/

echo ${xsec}

line2=${line2}
fout=${fout}
isdata=${isdata}
nowe=${nowe}
xsec=${xsec}
dout=${dout}
dout_f=${dout_f}
sample=${sample}
lumi=${lumi}
dataset=${dataset}
version=${version}

#Run the job
echo "Executing .././NtupleAnalyzer --file ${line2} --outfile ${dout_f}/${fout} --isdata ${isdata} --doSystCombine ${doSystCombine} --nowe ${nowe} --xsec ${xsec} --lumi ${lumi} --nmax ${nmax} --dataset ${dataset}"
${dout}/../Analyzer --file ${line2} --outfile ${dout_f}/${fout} --isdata ${isdata} --doSystCombine ${doSystCombine} --nowe ${nowe} --xsec ${xsec} --lumi ${lumi} --nmax ${nmax} --tree Nt --dataset ${dataset}

#Move output file to from scratch1 to grid, via script (args : output filename, and destination)
echo "Will move file [${dout_f}/${fout}.root] to [/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/${version}/${dataset}/.]"
${dout}/./Move_File_ToGrid.sh ${dout_f}/${fout}.root ${dataset} ${version}
