#!/bin/bash

#root -l -b -q "ProcessTree.C()" #produces histograms from tree
#root -l -b -q "ComputeV2.C()" #from 3D histograms produces 1D histograms of mass and V2 
for i in `seq 0 10` #OO centrality
    do
    root -l -b -q "FitV2orPol.C(1,0,1,0,$i)"
    done