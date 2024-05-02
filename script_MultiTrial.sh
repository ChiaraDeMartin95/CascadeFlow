#!/bin/bash

#for j in 0
for j in `seq 1 8` #number of trials
do
    root -l -b -q "ProcessTree.C($j)" #produces histograms from tree
    root -l -b -q "ComputeV2.C($j)" #from 3D histograms produces 1D histograms of mass and V2 

    for i in `seq 0 7` #centrality bins
    do
    root -l -b -q "FitV2.C($j, $i)" #fit of inv mass and V2 histograms    
    done
done