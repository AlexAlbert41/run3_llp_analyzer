# Instructions for C++ Analyzer, condor, and golden JSON

Follow the Setup instructions in the main README on the cmslpc. The OS doesn't matter, since the code will be run on an EL7 image instead.

---

## Running Analyzer Locally

The code will need to be compiled and run in an el7 image on the CMSLPC. See the instructions in the plotting script readme for opening the image.
Run `make`. If you get errors about something called "Fastjet", let me know :)
To create a single output ROOT file from a single MDSNANO input, you can run a command in the image like:
`./bin/Runllp_MuonSystem_CA_TrigEff_mdsnano Merged_Cache_InputLists/Muon0-Run2024B-PromptReco-v1-AOD.txt --isData -f=restestTrigEffMDSNANO.root -l=Summer24`
The `Merged_Cache_InputLists` directory houses the txt files that list the input MDSNANO files by era. I recommend making your own txt files with just a few input files that you can run over locally - you don't want to run over a whole era without condor.

## Running Over Condor

To efficiently run over all of the input files, you can use HTCondor, which can create many jobs for each Muon/Era that will run in parallel on the LPC condor nodes. An example command to run on condor is:
`python3 scripts_condor/submit_condor_LPC_TnP.py Muon0-Run2024B-PromptReco-v1-AOD TrigEff_mdsnano 2024_retest 2024_retest`
You can check the script to see exactly what the different command line arguments do. The output path is partially hardcoded to my directory within the shared eos area, so I would fix that before running. One subtlety is that condor is only available on el8, so even though `make` should be ran in el7 along with the local tests, the condor command needs to happen in el8. The os of the condor machines is specified in the script, so if we upgrade everything to work with a newer os, this will need to be modified.

## Hadd and Golden JSON

After all of your condor jobs are done, you can hadd the individual output files. I do this with a different condor script that I use for my current projects. When I was working on the trigger, I just did it by hand. If by hand proves annoying, let me know and I will share the code I use now, but may have to be tweaked for this data.

Lastly, we run the hadd'ed data over the golden json, which filters out events that happened during "bad" runs. Take a look at this repo, particularly this script that needs to be run.
https://github.com/AlexAlbert41/RazorCommon/blob/master/Tools/bin/run_goldenJSONLPC.sh
This repo needs to be installed in the <br> src <br> of your CMSSW directory, NOT within run3_llp_analyzer. Otherwise, things break. You will have to change the input/output directories and so forth. I also used this code most recently with a newer os/CMSSW, so I can't promise compatibility with CMSSW10. 
