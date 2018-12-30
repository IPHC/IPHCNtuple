#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

./../Analyzer \
--file input_synchro.txt \
--tree Nt \
--outfile ./output_tHq_MC \
--nmax -1  \
--isdata 0 \
--doSystCombine 0 \
--nowe 1 \
--xsec 1 \
--lumi 1 \
--dataset synchro



#THQ
#xsec = 0.2568
#SWE = 3361686
#--dataset ZZTo2L2Nu_13TeV_powheg_pythia8 #NB : if different from thq, ... => Training sel not filled
