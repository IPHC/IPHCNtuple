#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

./../Analyzer \
--file input.txt \
--tree Nt \
--outfile ./output_tHq_MC \
--nmax -1  \
--isdata 0 \
--doSystCombine 0 \
--nowe 6125476 \
--xsec 0.2568 \
--lumi 41500 \
--dataset THQ_ctcvcp_4f_Hincl_13TeV_madgraph_pythia8


#====================================
#====================================
#THQ_ctcvcp
#xsec = 0.2568
#SWE = 6125476
