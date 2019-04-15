#!/bin/env zsh

cdir=$(pwd)/..
export LD_LIBRARY_PATH=${cdir}:${cdir}/obj:$LD_LIBRARY_PATH

./NtupleProducer \
--file input.txt \
--outfile output \
--tree FlatTree/tree \
--nmax 20000 \
--isdata 0 \
--sync 0
