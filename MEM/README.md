#INSTALLATION

#See https://twiki.cern.ch/twiki/bin/viewauth/CMS/IPHCMEMCPP

#LHAPDF6 needs to be installed and LHAPDF environmental variables set.

cd CMSSW_8_0_23/src

cmsenv

git clone -b ttZDelphes https://github.com/IPHC/IPHCNtuple.git

cd IPHCNtuple/MEM/Madgraph

./bash SetupMadgraph.sh

cd ../Minimizer

root -l createLibSubGradient.C

make

cd ../src

make

cd ../test

#Need to update MEMEXECDIR

#Need to update config.cfg 

./test root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v3/output_Delphes_TTZselWithGen_4000k.root 0 1 \
 --MEMRun JobsReco/config.cfg

