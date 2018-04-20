#!/bin/bash

eosdir="/store/user/chanon/TZQ/TestNtuplesV7_Syst"

cmsLs ${eosdir} > l

for i in `cat l`
do

echo root://eoscms//eos/cms${eosdir}/$i

done

