#!/bin/sh

export X509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/
export SCRAM_ARCH=slc6_amd64_gcc630
eval `scramv1 runtime -sh`
cd -

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${dout}/../

echo "Executing ./Make_PU_Histograms.exe ..."

que="cms_local_mdm" #reserved slots, faster, <4h jobs only

rootdir="/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/"

logdir="${rootdir}log_PU"
mkdir ${logdir}

exepath="${rootdir}Make_PU_Histograms.exe"

${exepath}
