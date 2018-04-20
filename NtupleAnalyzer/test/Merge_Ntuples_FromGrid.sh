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

  if [[ $1 != "" && $line != *$1* ]]; then
	continue
  fi

  d1=$(echo $line)
  liMC2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}') #lists all FILENAMES in dir
  echo $liMC2 | while read line2
  do
        f1=$(echo $line2)
        file=$(echo ${fpathMCXRD}${d1}/${f1})
        echo "${file}" >> /tmp/tempMC.txt
  done

  echo "hadd -f $outDir/merged_$line.root " > /tmp/tmp_ScriptMerging_$line.sh
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

  if [[ $1 != "" && $line != *$1* ]]; then
	continue
  fi

  echo $line
  d1=$(echo $line)
  liDATA2=$(/usr/bin/rfdir ${fpath}${d1} | awk '{print $9}')
  echo $liDATA2 | while read line2
  do
        f1=$(echo $line2)
        file=$(echo ${fpathDATAXRD}${d1}/${f1})
        echo "${file}" >> /tmp/tempDATA.txt
  done

  echo "hadd -f $outDir/merged_$line.root " > /tmp/tmp_ScriptMerging_$line.sh
  cat /tmp/tempDATA.txt >> /tmp/tmp_ScriptMerging_$line.sh
  tr '\n' ' ' < /tmp/tmp_ScriptMerging_$line.sh > /tmp/ScriptMerging_$line.sh
  rm -f /tmp/tmp_ScriptMerging_$line.sh
  chmod 755 /tmp/ScriptMerging_$line.sh
  /tmp/ScriptMerging_$line.sh
  rm -f /tmp/ScriptMerging_$line.sh


  rm -f /tmp/tempDATA.txt
done


if [[ $1 != "" ]]; then
  exit 1
fi


#--- FURTHER HADDS & RENAME NTUPLES ---
echo ''
echo ''
echo ''
echo '=== FURTHER MERGING AND RENAMING ==='

# hadd -f $outDir/DY.root $outDir/merged_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root $outDir/merged_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root
# rm $outDir/merged_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root
# rm $outDir/merged_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root

hadd -f $outDir/DY.root $outDir/merged_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root $outDir/merged_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root

rm $outDir/merged_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f $outDir/TTbar_DiLep.root $outDir/merged_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root $outDir/merged_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root
rm $outDir/merged_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root
rm $outDir/merged_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root


hadd -f $outDir/TTbar_SingleLep.root $outDir/merged_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root $outDir/merged_TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root $outDir/merged_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root $outDir/merged_TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root
rm $outDir/merged_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root
rm $outDir/merged_TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root
rm $outDir/merged_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root
rm $outDir/merged_TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root


hadd -f $outDir/ttZ.root $outDir/merged_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root $outDir/merged_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2_v1_MINIAODSIM.root $outDir/merged_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3_v1_MINIAODSIM.root $outDir/merged_TTLL_m1to10.root
rm $outDir/merged_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root
rm $outDir/merged_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2_v1_MINIAODSIM.root
rm $outDir/merged_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3_v1_MINIAODSIM.root
rm $outDir/merged_TTLL_m1to10.root


hadd -f $outDir/ttG.root $outDir/merged_TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root $outDir/merged_TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root
rm $outDir/merged_TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root
rm $outDir/merged_TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root


hadd -f $outDir/tG.root $outDir/merged_TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root $outDir/merged_TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root
rm $outDir/merged_TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root
rm $outDir/merged_TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root

hadd -f $outDir/ZG.root $outDir/merged_ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root $outDir/merged_ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root
rm $outDir/merged_ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1_MINIAODSIM.root
rm $outDir/merged_ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM.root

hadd -f $outDir/GammaConv.root $outDir/ttG.root $outDir/ZG.root $outDir/tG.root $outDir/merged_WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root
rm $outDir/ttG.root
rm $outDir/ZG.root
rm $outDir/tG.root
rm $outDir/merged_WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f $outDir/ST.root $outDir/merged_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root $outDir/merged_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root $outDir/merged_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root
rm $outDir/merged_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root
rm $outDir/merged_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root
rm $outDir/merged_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root


mv $outDir/merged_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root $outDir/WZ.root
mv $outDir/merged_ST_tWll_5f_LO_13TeV-MadGraph-pythia8.root $outDir/tWZ.root
mv $outDir/merged_THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1.root $outDir/tHq.root
mv $outDir/merged_THW_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1.root $outDir/tHW.root
mv $outDir/merged_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root $outDir/ttW.root
mv $outDir/merged_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root $outDir/WWZ.root
mv $outDir/merged_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root $outDir/WZZ.root
mv $outDir/merged_ZZTo4L_13TeV_powheg_pythia8.root $outDir/ZZ.root
mv $outDir/merged_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root $outDir/ZZZ.root
mv $outDir/merged_TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root $outDir/TTTT.root
mv $outDir/merged_tZq_ll_4f_13TeV-amcatnlo-herwigpp.root $outDir/tZq.root
mv $outDir/merged_ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix.root $outDir/ttH.root
mv $outDir/merged_ttWJets_13TeV_madgraphMLM.root $outDir/ttW_MadgraphMLM.root
mv $outDir/merged_ttZJets_13TeV_madgraphMLM.root $outDir/ttZ_MadgraphMLM.root
mv $outDir/merged_WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root $outDir/WW_EWK.root
mv $outDir/merged_WW_DoubleScattering_13TeV-pythia8.root $outDir/WW_DPS.root
mv $outDir/merged_WWTo2L2Nu_13TeV-powheg.root $outDir/WWTo2L2Nu.root
mv $outDir/merged_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root $outDir/WJets.root
# mv $outDir/merged_TTLL_m1to10.root $outDir/TTLL_m1to10.root
# mv $outDir/merged_WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root $outDir/WG.root



hadd -f $outDir/SingleElectron.root $outDir/*SingleElectron*
hadd -f $outDir/SingleMuon.root $outDir/*SingleMuon*
hadd -f $outDir/MuonEG.root $outDir/*MuonEG*
hadd -f $outDir/DoubleEG.root $outDir/*DoubleEG*
hadd -f $outDir/DoubleMuon.root $outDir/*DoubleMuon*
hadd -f $outDir/DATA.root $outDir/DoubleMuon.root $outDir/DoubleEG.root $outDir/MuonEG.root $outDir/SingleElectron.root $outDir/SingleMuon.root
rm $outDir/*SingleElectron*2016*.root
rm $outDir/*SingleMuon*2016*.root
rm $outDir/*DoubleEG*2016*.root
rm $outDir/*DoubleMuon*2016*.root
rm $outDir/*MuonEG*2016*.root
