#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict

sample = sys.argv[1]
analyzer_str = sys.argv[2]
output_end = sys.argv[3]
log_name = sys.argv[4]

log_dir = "log_{}".format(log_name)
submit_dir = "submit_{}".format(log_name)
os.system("mkdir -p {}".format(log_dir))
os.system("mkdir -p {}".format(submit_dir))
executable = "scripts_condor/runAnalyzer.sh"
#analyzer = 'llp_MuonSystem_CA_TnP'
filesPerJob = 10
ntupler_version = 'V1p19/Data2023/'

#ntupler_version = "V1p19/MC_Summer22EE/v1/sixie/"


if analyzer_str=="TnP":analyzer="llp_MuonSystem_CA_TnP"
elif analyzer_str=="TnP_mdsnano":analyzer="llp_MuonSystem_CA_TnP_mdsnano"
elif analyzer_str=="TnP_noClusters":analyzer="llp_MuonSystem_CA_TnP_noClusters"
elif analyzer_str=="TrigEff":analyzer="llp_MuonSystem_CA_TrigEff"
elif analyzer_str=="TrigEff_mdsnano":analyzer="llp_MuonSystem_CA_TrigEff_mdsnano"
else:print("analyzer string not recognized, exiting ...");exit()
analyzer_version = 'v11'
#outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/{0}/{1}/".format(ntupler_version, analyzer_version)
#outputDirectoryBase="/store/group/lpclonglived/amalbert/Data_MC_Comp_TnP/results_from_cache_noSkim/{}/".format(output_end)
outputDirectoryBase="/store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/{}/".format(output_end)
HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')
#Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"
Analyzer_DIR = ''
#datasetListDir = Analyzer_DIR + "lists/displacedJetMuonNtuple/{}/".format(ntupler_version)
datasetListDir = Analyzer_DIR + "Merged_Cache_InputLists/"


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
print(sample)
if "Run2022E" in sample or "Run2022F" in sample or "Run2022G" in sample or "MC_Summer22EE" in sample:
    analyzerTag = "Summer22EE"
    jetVetoMap = "Summer22EE_23Sep2023_RunEFG_v1.root"
    pileupWeights = "PileupReweight_Summer22EE.root"
elif "Run2022C" in sample or "Run2022D" in sample or "MC_Summer22" in sample:
    analyzerTag = "Summer22"
    jetVetoMap = "Summer22_23Sep2023_RunCD_v1.root"
    pileupWeights = "PileupReweight_Summer22.root"
elif "Run2023D" in sample or "MC_Summer23BPix" in sample:
    analyzerTag = "Summer23BPix"
    jetVetoMap = "Summer23BPixPrompt23_RunD_v1.root"
    pileupWeights = "PileupReweight_Summer23BPix.root"
elif "Run2023B" in sample or "Run2023C" in sample or "MC_Summer23" in sample:
    analyzerTag = "Summer23"
    jetVetoMap = "Summer23Prompt23_RunC_v1.root"
    pileupWeights = "PileupReweight_Summer23.root"
elif "Run2024" in sample or "MC_Summer24" in sample:
    analyzerTag = "Summer24"
    jetVetoMap = "Winter24Prompt24_2024BCDEFGHI.root"
    pileupWeights = "PileupReweight_Summer24.root"
else:
    print("Couldn't find a valid analyzer tag. Exiting ...")
    exit()


#year = datasetList[sample][0]
#isData = datasetList[sample][1]
if "MC" in sample:isData="no"
else:isData="--isData"
#####################################
#Create Condor JDL file
#####################################
#os.system("rm -f submit/{}_{}_Job*.jdl".format(analyzer, sample))
#os.system("rm -f log/{}_{}_Job*".format(analyzer, sample))


jdl_file="{}/{}_{}_{}.jdl".format(submit_dir, analyzer, sample.replace("/","_"), maxjob)

tmpCondorJDLFile = open(jdl_file,"w")

tmpCondorJDLFile.write("Universe = vanilla \n")
tmpCondorJDLFile.write("Executable = {} \n".format(executable))
print("Arguments = {} {} {} {} $(ProcId) {} {} {} {} {}/ \n"\
                        .format(analyzer, inputfilelist, isData, filesPerJob, maxjob, outputDirectory, analyzerTag, CMSSW_BASE, HOME))
tmpCondorJDLFile.write("Arguments = {} {} {} {} $(ProcId) {} {} {} {} {}/ \n"\
                        .format(analyzer, inputfilelist, isData, filesPerJob, maxjob, outputDirectory, analyzerTag, CMSSW_BASE, HOME))

tmpCondorJDLFile.write("Log = {}/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).log \n".format(log_dir,analyzer, sample.replace("/","_"), maxjob))
tmpCondorJDLFile.write("Output = {}/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).out \n".format(log_dir, analyzer, sample.replace("/","_"), maxjob))
tmpCondorJDLFile.write("Error = {}/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).err \n".format(log_dir, analyzer, sample.replace("/","_"), maxjob))

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
tmpCondorJDLFile.write("Transfer_Input_Files = RazorRun,{},bin/Run{},DNNevaluation/EvaluateDNN.py,DNNevaluation/training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5, data/JetVetoMap/{}, data/PileupWeights/{}\n".format(inputfilelist,analyzer,jetVetoMap,pileupWeights))
tmpCondorJDLFile.write("Queue {} \n".format(maxjob))
#tmpCondorJDLFile.write("Queue {} \n".format(2))
tmpCondorJDLFile.close()

os.system("condor_submit {} --batch-name {}".format(jdl_file, sample))
