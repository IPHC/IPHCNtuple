#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

./../Analyzer \
--file input.txt \
--tree Nt \
--outfile ./output_tHq_MC \
--nmax 100  \
--isdata 0 \
--doSystCombine 0 \
--nowe 3495652 \
--xsec 0.7927 \
--lumi 35900 \

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
