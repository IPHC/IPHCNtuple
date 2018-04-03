#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

./../Analyzer \
--file input.txt \
--tree Nt \
--outfile ./output_Data \
--nmax 1000000  \
--isdata 1 \
--doSystCombine 0 \
--nowe -1 \
--xsec -1 \
--lumi 35900 \
