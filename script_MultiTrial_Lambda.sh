#!/bin/bash
#for j in `seq 0 199` #topo variations
for j in `seq 166 199` #topo variations
do
root -l -b -q "ComputeV2.C($j)" #from 3D histograms produces 1D histograms of mass and V2 
for i in `seq 0 10` #centrality
    do
     root -l -b -q "FitV2orPol.C(1, 0, 1, $j, $i)" #fit of inv mass and V2 histograms (j is for the mass cut)
    done
done

#for i in `seq 0 10` #centrality
#  do
#     root -l -b -q "MultiTrial.C($i, 3)" 
#  done
