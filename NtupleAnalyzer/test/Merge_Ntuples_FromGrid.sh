#!/bin/env zsh
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

# version="ttH2017"
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

#"$line" == *"madgraph_pythia8_Fall17"* ??

  if [[ $line == *"ZZTo2L2Nu"* || $line == *"WGToLNuG"* ]]; then
  echo "--> SKIPPED"
  continue
  fi

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

  #NEW -- CALL HERE A C++ SCRIPT TO GET THE HISTOGRAM CONTAINING THE SUMS OF WEIGHTS, AND STORE IT IN THE MERGED NTUPLE !
  ./Store_SumsOfWeights_inMergedNtuple.exe $line

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
#DY_LO
hadd -f $outDir/DY_LO.root $outDir/merged_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_v1_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_RECOSIMstep_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root

#DY_fxfx
hadd -f $outDir/DY_fxfx.root $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root

#DY_M4to50
hadd -f $outDir/DY_M4to50.root $outDir/merged_DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root $outDir/merged_DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root

rm $outDir/merged_DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root

#DYxJets
hadd -f $outDir/DYxJets_LO.root $outDir/merged_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root $outDir/merged_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root $outDir/merged_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root $outDir/merged_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root $outDir/merged_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM.root $outDir/merged_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v11_ext1_v1_MINIAODSIM.root $outDir/merged_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root
rm $outDir/merged_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v1_MINIAODSIM.root
rm $outDir/merged_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v10_v2_MINIAODSIM.root
rm $outDir/merged_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAOD_94X_mc2017_realistic_v11_ext1_v1_MINIAODSIM.root
rm $outDir/merged_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root

#FIXME -- merge DY_LO (=DYJets_M50) with DYxJets (=DYXXXJets_M50) ? #can not sum them ?
# hadd -f $outDir/DYxJets.root $outDir/DYxJets_LO.root $outDir/DY_LO.root
# rm $outDir/DYxJets_LO.root
# rm $outDir/DY_LO.root

hadd -f $outDir/ttZ.root $outDir/merged_TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/merged_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root
rm $outDir/merged_TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8.root
rm $outDir/merged_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root

hadd -f $outDir/TTJets_SemiLep_MLM.root $outDir/merged_TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8.root

# hadd -f $outDir/WJets.root $outDir/merged_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
# rm $outDir/merged_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
# rm $outDir/merged_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
# rm $outDir/merged_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root

# hadd -f $outDir/VHToNonbb.root $outDir/merged_VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_Fall17.root $outDir/merged_VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root
# rm $outDir/merged_VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root
# rm $outDir/merged_VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_Fall17.root

hadd -f $outDir/tWZ.root $outDir/merged_ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8.root $outDir/merged_ST_tWll_5f_LO_TuneCP5_PSweights_13TeV_madgraph_pythia8_Fall17.root
rm $outDir/merged_ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8.root
rm $outDir/merged_ST_tWll_5f_LO_TuneCP5_PSweights_13TeV_madgraph_pythia8_Fall17.root

hadd -f $outDir/ST_tWchan.root $outDir/merged_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root $outDir/merged_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root
rm $outDir/merged_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root
rm $outDir/merged_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root


hadd -f $outDir/WxJets.root $outDir/merged_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/merged_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root
rm $outDir/merged_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root

#//--------------------------------------------

#FCNC
mv $outDir/merged_ST_FCNC-TH_Tleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8.root $outDir/tH_ST_hut_FCNC.root
hadd -f $outDir/tH_TT_hut_FCNC.root $outDir/merged_TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8.root $outDir/merged_TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8.root
rm $outDir/merged_TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8.root
rm $outDir/merged_TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hut-MadGraph5-pythia8.root

mv $outDir/merged_ST_FCNC-TH_Tleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8.root $outDir/tH_ST_hct_FCNC.root
hadd -f $outDir/tH_TT_hct_FCNC.root $outDir/merged_TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8.root $outDir/merged_TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8.root
rm $outDir/merged_TT_FCNC-TtoHJ_aTleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8.root
rm $outDir/merged_TT_FCNC-aTtoHJ_Tleptonic_HToWWZZtautau_eta_hct-MadGraph5-pythia8.root



#//--------------------------------------------
mv $outDir/merged_THQ_4f_Hincl_13TeV_madgraph_pythia8.root $outDir/tHq_old.root
mv $outDir/merged_THW_5f_Hincl_13TeV_madgraph_pythia8.root $outDir/tHW_old.root
mv $outDir/merged_THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8.root $outDir/tHq.root
mv $outDir/merged_THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8.root $outDir/tHW.root

mv $outDir/merged_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root $outDir/ttG.root
mv $outDir/merged_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root $outDir/TTbar_DiLep.root
mv $outDir/merged_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root $outDir/TTbar_DiLep_PSweights.root
mv $outDir/merged_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root $outDir/TTbar_SemiLep.root
mv $outDir/merged_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root $outDir/TTbar_SemiLep_PSweights.root
mv $outDir/merged_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root $outDir/ttW.root
mv $outDir/merged_TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8.root $outDir/ttW_PSweights.root
mv $outDir/merged_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WWW.root
mv $outDir/merged_WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WWZ.root
mv $outDir/merged_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root $outDir/WZ.root
mv $outDir/merged_WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WZZ.root
mv $outDir/merged_ZZZ_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/ZZZ.root
mv $outDir/merged_tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8.root $outDir/tZq.root
mv $outDir/merged_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8.root $outDir/ttH.root
mv $outDir/merged_TTWW_TuneCP5_13TeV-madgraph-pythia8.root $outDir/ttWW.root
mv $outDir/merged_TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8.root $outDir/tGJets.root
mv $outDir/merged_TTTT_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/TTTT.root
mv $outDir/merged_ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8.root $outDir/ttH_LO_nonbb.root
mv $outDir/merged_ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8.root $outDir/ttW_LO.root
mv $outDir/merged_ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8.root $outDir/ttZ_LO.root
mv $outDir/merged_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root $outDir/TTJets.root
mv $outDir/merged_TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/TTJets_DiLep_MLM.root
mv $outDir/merged_ZZTo4L_13TeV_powheg_pythia8.root $outDir/ZZ.root
# mv $outDir/merged_ZZTo4L_13TeV_powheg_pythia8_RunIIFall17MiniAODv2_PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_v1_MINIAODSIM.root $outDir/ZZ.root
mv $outDir/WZG_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WZG.root
mv $outDir/merged_WpWpJJ_EWK-QCD_TuneCP5_13TeV-madgraph-pythia8.root $outDir/WpWp.root
mv $outDir/merged_WW_DoubleScattering_13TeV-pythia8_TuneCP5.root $outDir/WW_DPS.root
mv $outDir/merged_TTWH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/ttWH.root
mv $outDir/merged_TTTW_TuneCP5_13TeV-madgraph-pythia8.root $outDir/tttW.root
mv $outDir/merged_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root $outDir/GGHZZ4L.root
mv $outDir/merged_WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8.root $outDir/WG.root
mv $outDir/merged_WZG_TuneCP5_13TeV-amcatnlo-pythia8.root $outDir/WZG.root
mv $outDir/merged_ttH_M125_TuneCP5_13TeV-powheg-pythia8.root $outDir/ttH_LO.root
mv $outDir/merged_TprimeBToTH_M-600_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M600.root
mv $outDir/merged_TprimeBToTH_M-650_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M650.root
mv $outDir/merged_TprimeBToTH_M-700_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M700.root
mv $outDir/merged_TprimeBToTH_M-800_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M800.root
mv $outDir/merged_TprimeBToTH_M-900_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M900.root
mv $outDir/merged_TprimeBToTH_M-1000_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M1000.root
mv $outDir/merged_TprimeBToTH_M-1100_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M1100.root
mv $outDir/merged_TprimeBToTH_M-1200_LH_TuneCP5_13TeV-madgraph-pythia8.root $outDir/VLQ_M1200.root
mv $outDir/merged_VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root $outDir/VHToNonbb.root




#//--------------------------------------------
hadd -f $outDir/DATA.root $outDir/*SingleElectron* $outDir/*SingleMuon* $outDir/*MuonEG* $outDir/*DoubleEG* $outDir/*DoubleMuon*

rm $outDir/*SingleElectron*2017*.root
rm $outDir/*SingleMuon*2017*.root
rm $outDir/*DoubleEG*2017*.root
rm $outDir/*DoubleMuon*2017*.root
rm $outDir/*MuonEG*2017*.root
# //--------------------------------------------
# hadd -f $outDir/SingleElectron.root $outDir/*SingleElectron*
# hadd -f $outDir/SingleMuon.root $outDir/*SingleMuon*
# hadd -f $outDir/MuonEG.root $outDir/*MuonEG*
# hadd -f $outDir/DoubleEG.root $outDir/*DoubleEG*
# hadd -f $outDir/DoubleMuon.root $outDir/*DoubleMuon*

# hadd -f $outDir/DATA.root $outDir/DoubleMuon.root $outDir/DoubleEG.root $outDir/MuonEG.root $outDir/SingleElectron.root $outDir/SingleMuon.root

# rm $outDir/*SingleElectron*2017*.root
# rm $outDir/*SingleMuon*2017*.root
# rm $outDir/*DoubleEG*2017*.root
# rm $outDir/*DoubleMuon*2017*.root
# rm $outDir/*MuonEG*2017*.root

# rm $outDir/SingleElectron.root
# rm $outDir/SingleMuon.root
# rm $outDir/DoubleEG.root
# rm $outDir/DoubleMuon.root
# rm $outDir/MuonEG.root
#//--------------------------------------------
