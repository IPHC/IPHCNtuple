#!/bin/env zsh

cdir=$(pwd)/../
export LD_LIBRARY_PATH=${cdir}:$LD_LIBRARY_PATH

infl="input.txt"

./NtupleProducer      \
--file ${infl}        \
--outfile output \
--tree FlatTree/tree  \
--nmax 20000              \
--isdata 1
