#!/usr/bin/python

import numpy as np
import pandas as pd
import uproot
#import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
#import mplhep as hep
#import pickle
#import glob
import ROOT as rt
#import coffea
import awkward as ak
#from coffea import hist, processor
#from coffea.nanoevents.methods import candidate
#from coffea.nanoevents.methods import vector

#sys.path.append("/uscms/home/amalbert/nobackup/CMSSW_14_1_0_pre4/src/RazorCommon/Tools/bin")
#import importlib
#import getMuonScaleFactor

input_MC_file_23 = "/uscms/home/amalbert/nobackup/CMSSW_10_6_20/src/run3_llp_analyzer/Merged_Cache_InputLists/MC_Summer23/DYto2Mu_MLL-50to120.txt"
input_MC_file_23BPix = "/uscms/home/amalbert/nobackup/CMSSW_10_6_20/src/run3_llp_analyzer/Merged_Cache_InputLists/MC_Summer23BPix/DYto2Mu_MLL-50to120.txt"
input_MC_file_22 = "/uscms/home/amalbert/nobackup/CMSSW_10_6_20/src/run3_llp_analyzer/Merged_Cache_InputLists/MC_Summer22/DYto2Mu_MLL-50to120_keepMDSHits.txt"
input_MC_file_22EE = "/uscms/home/amalbert/nobackup/CMSSW_10_6_20/src/run3_llp_analyzer/Merged_Cache_InputLists/MC_Summer22EE/DYto2Mu_MLL-50to120_keepMDSHits.txt"
input_MC_file_24 = "/uscms/home/amalbert/nobackup/CMSSW_10_6_20/src/run3_llp_analyzer/Merged_Cache_InputLists/MC_Summer24/DYto2Mu_MLL-50to120.txt"

def computeWeightSum(input_MC_file):
    inputFileList = []
    with open(input_MC_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "Job217_of_226" in line: continue
            inputFileList.append(line)
    print("number of input files: ", len(inputFileList))
    total_files = len(inputFileList)
    increases = []
    cumulative_minus20 = 0
    weight_sum = 0
    for i in range(total_files):
        #if i<1080: continue
        if i%20==0:
            print("On file {} of total_files".format(i))
            print("sum of weights = ", weight_sum)
            print("weight increase from last 20: ",weight_sum-cumulative_minus20)
            increases.append(weight_sum-cumulative_minus20)
            cumulative_minus20 = weight_sum
        #print(inputFileList[i])
        #print(badFiles[0])
        #print(type(inputFileList[i]))
        #if inputFileList[i] in badFiles:
        #    continue
        myfile = inputFileList[i]
        #print(myfile)
        myfile = myfile[:-1]
        #print(myfile)
        myfile = uproot.open(myfile)
        f = myfile["Events"]
        #f = uproot.lazy(test_file_list[0]+":Events;1", "genWeight", "100 MB")
        arr = f["genWeight"].array(library="np")
        weight_sum+=np.sum(arr)

    return weight_sum

if __name__ == "__main__":
    print("Starting sum of weights")
    MC_Summer24_WeightSum = computeWeightSum(input_MC_file_24)
    print("MC Weight Sum for MC_Summer24 = ", MC_Summer24_WeightSum)
