#!/bin/bash
#for j in `seq 21 99` 
#do 
#root -l -b -q "ComputeV2.C($j)"
for i in `seq 0 9` #centrality bins
#for i in `seq 3 3`
 #do
 #for j in `seq 0 19` #BDT
do
  #for l in `seq 0 8` #mass
    #do
    #root -l -b -q "ProcessTHN.C($j)" #produces histograms from tree
    #root -l -b -q "ComputeV2.C($j)" #from 3D histograms produces 1D histograms of mass and V2 
    #root -l -b -q "ComputeEff.C(0, $j)"    
    root -l -b -q "FitV2orPol.C(1, 0, 1, 0, $i)" #fit of inv mass and V2 histograms (j is for the mass cut)
done
#done