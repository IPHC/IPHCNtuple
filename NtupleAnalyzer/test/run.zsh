#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

./../Analyzer \
--file input.txt \
--tree Nt \
--outfile ./output_tHq_MC \
--nmax 10000  \
--isdata 0 \
--doSystCombine 0 \
--lumi 41500 \
--nowe 9865010 \
--xsec 0.7927 \
--dataset THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8 \

#====================================
#====================================
#THQ_ctcvcp
#xsec = 0.7927
#nowe = 9865010

#TTWJets
#xsec = 0.2043
#nowe = 2678775

#ST_tWll
#xsec = 0.01123
#nowe = 927720
