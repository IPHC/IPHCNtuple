#!/bin/env zsh
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

# version="newlepMVA"
version="tHq2017"

fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/"$version"/"

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

#//--------------------------------------------
hadd -f $outDir/DY.root $outDir/merged_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root
rm $outDir/merged_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root


hadd -f $outDir/ttZ.root $outDir/merged_TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/merged_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root
rm $outDir/merged_TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8.root
rm $outDir/merged_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root


hadd -f $outDir/ST_tchan.root $outDir/merged_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root  $outDir/merged_ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root
rm $outDir/merged_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root
rm $outDir/merged_ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root

hadd -f $outDir/ST_tW.root $outDir/merged_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root  $outDir/merged_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root
rm $outDir/merged_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root
rm $outDir/merged_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root
#//--------------------------------------------


#//--------------------------------------------
mv $outDir/merged_THQ_4f_Hincl_13TeV_madgraph_pythia8.root $outDir/tHq.root
mv $outDir/merged_THQ_4f_Hincl_13TeV_madgraph_pythia8_Fall17.root $outDir/tHq_PrivateProd.root
mv $outDir/merged_THW_5f_Hincl_13TeV_madgraph_pythia8_Fall17.root $outDir/tHW_PrivateProd.root
mv $outDir/merged_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root $outDir/ttG.root
mv $outDir/merged_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root $outDir/TTbar_DiLep_PSweights.root
mv $outDir/merged_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root $outDir/TTbar_DiLep.root
mv $outDir/merged_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root $outDir/TTbar_Hadronic.root
mv $outDir/merged_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root $outDir/TTbar_Hadronic_PSweights.root
mv $outDir/merged_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root $outDir/TTbar_SemiLep.root
mv $outDir/merged_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root $outDir/TTbar_SemiLep_PSweights.root
mv $outDir/merged_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root $outDir/ttW.root
mv $outDir/merged_TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8.root $outDir/ttW_PSweights.root
mv $outDir/merged_WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8.root $outDir/WW.root
mv $outDir/merged_WWTo2L2Nu_DoubleScattering_13TeV-herwigpp.root $outDir/WW_DS.root
mv $outDir/merged_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WWW.root
mv $outDir/merged_WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WWZ.root
mv $outDir/merged_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root $outDir/WZ.root
mv $outDir/merged_WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WZZ.root
mv $outDir/merged_ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root $outDir/ZZ_otherVersion.root
mv $outDir/merged_ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM.root $outDir/ZZ.root
mv $outDir/merged_ZZZ_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/ZZZ.root
mv $outDir/merged_tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8.root $outDir/tZq.root
mv $outDir/merged_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8.root $outDir/ttH.root
mv $outDir/merged_ST_tWll_5f_LO_TuneCP5_PSweights_13TeV_madgraph_pythia8_Fall17.root $outDir/tWZ.root
#//--------------------------------------------

#//--------------------------------------------
hadd -f $outDir/SingleElectron.root $outDir/*SingleElectron*
hadd -f $outDir/SingleMuon.root $outDir/*SingleMuon*
hadd -f $outDir/MuonEG.root $outDir/*MuonEG*
hadd -f $outDir/DoubleEG.root $outDir/*DoubleEG*
hadd -f $outDir/DoubleMuon.root $outDir/*DoubleMuon*
hadd -f $outDir/DATA.root $outDir/DoubleMuon.root $outDir/DoubleEG.root $outDir/MuonEG.root $outDir/SingleElectron.root $outDir/SingleMuon.root
rm $outDir/*SingleElectron*2017*.root
rm $outDir/*SingleMuon*2017*.root
rm $outDir/*DoubleEG*2017*.root
rm $outDir/*DoubleMuon*2017*.root
rm $outDir/*MuonEG*2017*.root
#//--------------------------------------------
