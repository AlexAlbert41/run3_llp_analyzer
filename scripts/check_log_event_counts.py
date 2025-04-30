import sys
import os
import glob
import numpy as np

sample = sys.argv[2]
directory = sys.argv[1]

directory = "/uscms/home/amalbert/nobackup/CMSSW_10_6_20/src/run3_llp_analyzer/"+directory

pathList = glob.glob(f'{directory}/*{sample}*.out')
print(f"iterating through {len(pathList)} files")

sumEvents = 0
finalEvents = 0
eventsTwoTagged = 0
efficiencies = []
twotagefficiencies = []
for num, file in enumerate(pathList):
    if num%100==0:print(f"on file {num} of sample {sample}")
    f = open(file, 'r')
    lines = f.readlines()
    foundProcessed=False
    for line in lines:
        if "Processed" in line:
            foundProcessed=True
            newLine = line.split(" ")
            sumEvent = int(newLine[1])
            sumEvents+=sumEvent
        if "ZMuMu events selected" in line:
            newLine = line.split("=")
            finalEvent=int(newLine[1])
            finalEvents+=finalEvent
        if "Events Two Tagged" in line:
            newLine = line.split(":")
            eventTwoTagged=int(newLine[1])
            eventsTwoTagged+=eventTwoTagged
    if not foundProcessed:
        print(f"Data not processed for file {file}")
    efficiencies.append(finalEvent/sumEvent)
    twotagefficiencies.append(eventTwoTagged/sumEvent)
efficiencies_arr = np.array(efficiencies)
twotagefficiencies_arr = np.array(twotagefficiencies)
print(f"Average efficiency = {np.average(efficiencies_arr)*100}%")
print(f"Max efficiency = {np.max(efficiencies_arr)*100}%")
print(f"Min efficiency = {np.min(efficiencies_arr)*100}%")
print(f"Efficiency StDev = {np.std(efficiencies_arr)*100}%")
print(f"Average Fraction of events with two tags = {np.average(twotagefficiencies_arr)*100}%")
print(f"two tag inclusive events = {(finalEvents+eventsTwoTagged)}")
print(f"total two tag inclusive efficiencuy = {(finalEvents+eventsTwoTagged)/(2*(sumEvents))*100}%")
print(f"Processed {sumEvents} for sample {sample}")
print(f"{finalEvents} in output for sample {sample}")
        
