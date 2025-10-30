#!/bin/bash

#root -l -b -q "ProcessTree.C()" #produces histograms from tree
#root -l -b -q "ComputeV2.C()" #from 3D histograms produces 1D histograms of mass and V2 
#MultType = 1 for FT0M, isMB = 0, mul = i
for i in `seq 0 8`
    do
    #root -l -b -q "FitV2orPol.C(1,1,0,$i)"
    #root -l -b -q "FitV2orPol.C(1,0,0,$i)"
    root -l -b -q "FitV2orPol.C(1,1,0,$i)"
    #root -l -b -q "FitV2orPol.C(0,1,0,$i,2)" #XiMinus
    #root -l -b -q "FitV2orPol.C(0,0,0,$i,2)" #XiMinus
    #root -l -b -q "FitV2orPol.C(1,1,0,$i,3)" #XiPlus
    #root -l -b -q "FitV2orPol.C(1,0,0,$i,3)" #XiPlus
    #root -l -b -q "FitV2orPol.C(1,0,0,$i,3)" #XiPlus
    #root -l -b -q "FitV2orPol.C(1,1,0,$i,2)" #XiMinus
    #root -l -b -q "FitV2orPol.C(1,1,0,$i,3)" #XiPlus
    #root -l -b -q "FitV2orPol.C(0,0,0,$i,2)" #XiMinus
    #root -l -b -q "FitV2orPol.C(0,0,0,$i,3)" #XiPlus
    #root -l -b -q "FitV2orPol.C(0,1,0,$i,2)" #XiMinus
    #root -l -b -q "FitV2orPol.C(0,1,0,$i,3)" #XiPlus
    #root -l -b -q "FitV2orPol.C(1,0,0,$i,0)" #Xi v2
    done

#pt analysis, polfromLambda