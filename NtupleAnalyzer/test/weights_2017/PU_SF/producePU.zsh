#!/bin/env zsh

#pp inelastic xsec at 13 TeV // +- 5% variations
xsec=("67735" "71300" "74865")
tag=("Down" "Nom" "Up")

idx=1
for i in ${xsec}
do
echo $i
pileupCalc.py -i /home-pbs/ntonon/tHq/FlatTreeProducer_2017/CMSSW_9_4_3/src/IPHCFlatTree/FlatTreeProducer/test/PROD/GRL/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt \
--inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt \
--calcMode true \
--minBiasXsec ${i} \
--maxPileupBin 99 \
--numPileupBins 99 \
Pileup${tag[${idx}]}.root
idx=$[$idx+1]
done


#--maxPileupBin 50 \
#--numPileupBins 50 \
