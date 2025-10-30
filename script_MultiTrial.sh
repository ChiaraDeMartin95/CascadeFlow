#!/bin/bash
#for i in `seq 4 8` #centrality bins
#for i in `seq 3 3`
 #do
 for j in `seq 0 19` #BDT
  do
  #for l in `seq 0 8` #mass
    #do
    root -l -b -q "ProcessTHN.C($j)" #produces histograms from tree
    root -l -b -q "ComputeV2.C($j)" #from 3D histograms produces 1D histograms of mass and V2 
    #root -l -b -q "ComputeEff.C(0, $j)"    
    #root -l -b -q "FitV2orPol.C(1, 0, 0, $i, $j)" #fit of inv mass and V2 histograms (j is for the mass cut)
    #root -l -b -q "FitV2orPol.C(1, 0, 0, $i)" #fit of inv mass and V2 histograms (j is for the mass cut)
    #root -l -b -q "FitV2orPol.C(1, 1, $j, $i, $l)" #fit of inv mass and V2 histograms (j is for the mass cut)   
    #root -l -b -q "FitV2orPol.C(1, 0, 0, $i, 0, 1, 6)" #for Lambdas
    #root -l -b -q "MultiTrial.C($i, 2)" #fit of inv mass and V2 histograms    
    #echo "hello"
    #for j in `seq 0 21` #number of trials
    #do
    #for l in `seq 0 6` #number of mass interval trials
    #do
    #root -l -b -q "FitV2orPol.C(1, 1, $j, $i, $l)" #fit of inv mass and V2 histograms      
    #done
   done
#done