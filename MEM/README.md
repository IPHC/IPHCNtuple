#INSTALLATION

#See https://twiki.cern.ch/twiki/bin/viewauth/CMS/IPHCMEMCPP

#LHAPDF6 needs to be installed and LHAPDF environmental variables set.

cd CMSSW_8_0_20/src

cmsenv

git clone -b tZqEPS2017 https://github.com/IPHC/IPHCNtuple.git

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

./test root://eoscms//eos/cms/store/user/chanon/TZQ/TestNtuplesV8bTight/FCNCNTuple_tZqmcNLO.root 0 1 \
 --MEMRun JobsReco/config.cfg

