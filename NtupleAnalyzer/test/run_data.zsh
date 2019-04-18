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
--lumi 1 \
--nowe 1 \
--xsec 1 \
--dataset SingleMuon_Run2017B_17Nov2017_v1_MINIAOD \
