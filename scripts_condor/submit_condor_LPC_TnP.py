#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict


os.system("mkdir -p submit_data2")
os.system("mkdir -p log_data2")
executable = "scripts_condor/runAnalyzer.sh"
analyzer = 'llp_MuonSystem_CA_TnP'
filesPerJob = 10
ntupler_version = 'V1p19/Data2023/'

#ntupler_version = "V1p19/MC_Summer22EE/v1/sixie/"

sample = sys.argv[1]

analyzer_version = 'v11'
#outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/{0}/{1}/".format(ntupler_version, analyzer_version)
outputDirectoryBase="/store/user/amalbert/MDSTriggerEff/Run2023_condor_TnP_Output_fixTagPt_103124/"
HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')
#Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"
Analyzer_DIR = ''
#datasetListDir = Analyzer_DIR + "lists/displacedJetMuonNtuple/{}/".format(ntupler_version)
datasetListDir = Analyzer_DIR + "newMerge_pathToData/"


'''
#### Run on all signal samples ####
## year/isData
datasetList = OrderedDict()
samples = os.listdir(datasetListDir)
for s in samples:
    if not "EXOCSCCluster" in s:continue 
    if 'Data' in ntupler_version: datasetList[s.replace('.txt', '')] = ["2018", "yes"]
    else: datasetList[s.replace('.txt', '')] = ["2018", "no"]
############
'''



#for sample in datasetList.keys():

print("Preparing analyzer workflow for dataset: " + sample + "\n")

inputfilelist  = datasetListDir + sample +'.txt'
if not os.path.exists(inputfilelist):
    print("listfile: " + inputfilelist + " does not exist. skipping.")
    #continue
cmd = ["awk", "END{print NR}",inputfilelist]
nfiles = int(subprocess.check_output(cmd).decode("utf-8"))
maxjob=int(nfiles/filesPerJob)+1
mod=int(nfiles%filesPerJob)
if mod == 0: maxjob -= 1
outputDirectory = outputDirectoryBase + sample + "/"
analyzerTag = "test"
#year = datasetList[sample][0]
#isData = datasetList[sample][1]
isData="no"
#####################################
#Create Condor JDL file
#####################################
#os.system("rm -f submit/{}_{}_Job*.jdl".format(analyzer, sample))
#os.system("rm -f log/{}_{}_Job*".format(analyzer, sample))


jdl_file="submit_data2/{}_{}_{}.jdl".format(analyzer, sample, maxjob)

tmpCondorJDLFile = open(jdl_file,"w")

tmpCondorJDLFile.write("Universe = vanilla \n")
tmpCondorJDLFile.write("Executable = {} \n".format(executable))
print("Arguments = {} {} {} {} $(ProcId) {} {} {} {} {}/ \n"\
                        .format(analyzer, inputfilelist, isData, filesPerJob, maxjob, outputDirectory, analyzerTag, CMSSW_BASE, HOME))
tmpCondorJDLFile.write("Arguments = {} {} {} {} $(ProcId) {} {} {} {} {}/ \n"\
                        .format(analyzer, inputfilelist, isData, filesPerJob, maxjob, outputDirectory, analyzerTag, CMSSW_BASE, HOME))

tmpCondorJDLFile.write("Log = log_data2/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).log \n".format(analyzer, sample, maxjob))
tmpCondorJDLFile.write("Output = log_data2/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).out \n".format(analyzer, sample, maxjob))
tmpCondorJDLFile.write("Error = log_data2/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).err \n".format(analyzer, sample, maxjob))

tmpCondorJDLFile.write("+JobQueue=\"Short\" \n")
tmpCondorJDLFile.write("RequestMemory = 4096 \n")
tmpCondorJDLFile.write("RequestCpus = 1 \n")
tmpCondorJDLFile.write("RequestDisk = 4 \n")

tmpCondorJDLFile.write("+RunAsOwner = True \n")
tmpCondorJDLFile.write("+InteractiveUser = true \n")
tmpCondorJDLFile.write("+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7\" \n")
tmpCondorJDLFile.write('+SingularityBindCVMFS = True \n')
tmpCondorJDLFile.write("run_as_owner = True \n")
#tmpCondorJDLFile.write("x509userproxy = {}/x509_proxy \n".format(HOME))
tmpCondorJDLFile.write("should_transfer_files = YES \n")
tmpCondorJDLFile.write("when_to_transfer_output = ON_EXIT_OR_EVICT \n")
tmpCondorJDLFile.write("Transfer_Input_Files = RazorRun,{},bin/Runllp_MuonSystem_CA_TnP,DNNevaluation/EvaluateDNN.py,DNNevaluation/training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5\n".format(inputfilelist))
tmpCondorJDLFile.write("Queue {} \n".format(maxjob))
#tmpCondorJDLFile.write("Queue {} \n".format(2))
tmpCondorJDLFile.close()

os.system("condor_submit {} --batch-name {}".format(jdl_file, sample))
