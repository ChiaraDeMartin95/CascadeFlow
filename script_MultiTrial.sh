#!/bin/bash
#for j in `seq 1 20` 
#do
#    root -l -b -q "ComputeV2.C($j)" #from 3D histograms produces 1D histograms of mass and V2 
#done

for j in `seq 81 100` 
do
root -l -b -q "FitV2orPol_CB.C(0, 0, 0, 1, $j, 10)"
for i in `seq 0 10`
do
  root -l -b -q "FitV2orPol_CB.C(8, 1, 0, 1, $j, $i)"
done
done