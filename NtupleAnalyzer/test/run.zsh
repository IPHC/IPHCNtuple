#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

./../Analyzer \
--file input.txt \
--tree Nt \
--outfile ./output_tHq_MC \
--nmax -1  \
--isdata 1 \
--doSystCombine 0 \
--nowe 3495652 \
--xsec 0.7927 \
--lumi 35.9 \
--dataset DoubleEG_Run2016B_23Sep2016_v3_MINIAOD\

#./../Analyzer \
#--file input.txt \
#--tree Nt \
#--outfile ./output_tHq_MC \
#--nmax -1 \
#--isdata 0 \
#--doSystCombine 0 \
#--nowe 1 \
#--xsec 1 \
#--lumi 1 \
