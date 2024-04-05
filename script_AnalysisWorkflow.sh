#!/bin/bash

root -l -b -q "ProcessTree.C()" #produces histograms from tree
root -l -b -q "ComputeV2.C()" #from 3D histograms produces 1D histograms of mass and V2 
#MultType = 1 for FT0M, isMB = 0, mul = i
for i in `seq 0 7`
    do
    root -l -b -q "FitV2.C($i)" #fit of inv mass and V2 histograms    
    done