#!/bin/env zsh

cdir=$(pwd)/..
export LD_LIBRARY_PATH=${cdir}:${cdir}/obj:$LD_LIBRARY_PATH

infl="input.txt"
nmax=-1

./NtupleProducer \
--file ${infl} \
--outfile output \
--tree FlatTree/tree \
--nmax ${nmax} \
--isdata 0 \
--sync 0
