#!/bin/sh

hostname
echo "Job started"
date
start_time=`date +%s`
analysisType=$1
inputfilelist=$2
isData=$3
filePerJob=$4
jobnumber=$5
maxjob=$6
sample=${inputfilelist##*/}
sample=${sample%.txt}
outputfile=${sample}_Job${jobnumber}_of_${maxjob}.root
outputDirectory=$7
analyzerTag=$8
CMSSW_BASE=$9
homeDir=${10}
currentDir=`pwd`
echo "currentDir: ${currentDir}"
#user=${homeDir#*/data/}
#user=${homeDir#*/storage/user/}
#runDir=${currentDir}/${analyzerTag}/
runDir=${currentDir}/CMSSW_10_6_20/src/

rm -rf ${runDir}
#mkdir -p ${runDir}
echo ${CMSSW_BASE}
echo homeDir: ${homeDir}

:'
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	#setup cmssw
	ls -la
	cd CMSSW_10_6_20/src/
	#workDir=`pwd`
	#echo "entering directory: ${workDir}"
	ulimit -c 0
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	eval `scram runtime -sh`
	cd ${currentDir}
	#echo `which root`

	#cd ${runDir}
	#echo "entering directory: ${runDir}"
	#echo "${CMSSW_BASE}/src/run3_llp_analyzer/RazorRun"
'
#copying LPC commands

echo ${CMSSW_BASE}
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	export CWD=${PWD}
	export PATH=${PATH}:/cvmfs/cms.cern.ch/common/
	export SCRAM_ARCH=slc7_amd64_gcc700
	echo "PATH: $PATH"
	echo "SCRAM_ARCH: $SCRAM_ARCH"
	scramv1 project CMSSW CMSSW_10_6_20
	#cmsrel CMSSW_10_6_20
	echo "Inside $currentDir:"
	ls -lah
	#cp RazorRun CMSSW_10_6_20/src/run3_llp_analyzer/
	echo "CMSSW_BASE ${CMSSW_BASE}"
	#cd ${CMSSW_BASE}
	#cd CMSSW_10_6_20/src
	ls -lah

#cp RazorRun CMSSW_10_6_20/src
#cp Runllp_MuonSystem_CA_TnP CMSSW_10_6_20/src
#cp EvaluateDNN.py CMSSW_10_6_20/src
#cp training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5
	if [ -f RazorRun ]
	then
		#cp $CMSSW_BASE/src/run3_llp_analyzer/RazorRun ./
		#cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/EvaluateDNN.py ./
        #        cp $CMSSW_BASE/src/run3_llp_analyzer/DNNevaluation/*.h5 .
		cp RazorRun ${runDir}
		cp Runllp_MuonSystem_CA_TnP ${runDir}
		cp EvaluateDNN.py ${runDir}
		cp training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5 ${runDir}
		cp ${sample}.txt ${runDir}

		
		#get grid proxy
		export X509_USER_PROXY=${currentDir}/x509up_u57571
		echo "${currentDir}/x509up_u57571"
		voms-proxy-info

		cd ${runDir}
		echo "entering directory: ${runDir}"
		#cmsenv
		eval `scram runtime -sh`

		#run the job
		# echo "cat ${inputfilelist} | awk \"NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})\" > inputfilelistForThisJob_${jobnumber}.txt"
		cat ${sample}.txt | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
                echo "Running on these input files:"
                cat inputfilelistForThisJob_${jobnumber}.txt
                echo "************************************"
	
		echo "start copying files"
		while read line; do   xrdcp -f ${line} .; done < "inputfilelistForThisJob_${jobnumber}.txt"
		rm inputfilelistForThisJob_${jobnumber}.txt
		ls -ltr *ntupler*.root
		ls -ltr
		echo "now sending ls output to new file"
		ls *ntupler*.root > inputfilelistForThisJob_${jobnumber}.txt
		
		

		echo ""
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
		echo ""
		echo " "; echo "Starting razor run job now"; echo " ";
		if [ ${analysisType} == "MakeMCPileupDistribution" ]
		then
			echo "./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -f=${outputfile}"
			./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -f=${outputfile}
		else
			echo ./Runllp_MuonSystem_CA_TnP inputfilelistForThisJob_${jobnumber}.txt -d=${isData}  -f=${outputfile} -l=${analyzerTag}
			./Runllp_MuonSystem_CA_TnP inputfilelistForThisJob_${jobnumber}.txt  --isData  -f=${outputfile}
		fi

		echo ${outputfile}
		echo ${outputDirectory}
		ls *root > output.txt
		echo "Output ROOT files: "
		cat output.txt
		##^_^##
		echo "RazorRun_T2 finished"
		date

		echo "start DNN evaluation"
		source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh
		python EvaluateDNN.py --in_file ${outputfile}
		
		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to ${outputDirectory}"
		eval `scram unsetenv -sh`
		#mkdir -p ${outputDirectory}
		echo "trying to copy file"
		xrdcp ${outputfile} root://cmseos.fnal.gov/${outputDirectory}/${outputfile}
		echo "copied file"
		:'
		if [ -f ${outputDirectory}/${outputfile} ]
		then
			echo ${outputfile} "copied"
		else
			echo ${outputfile} "not copied"
		fi
		'
	else
		echo echo "WWWWYYYY ============= failed to access file RazorRun_T2, job anandoned"
	fi

else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
#rm -rf ${runDir}
echo "Job finished"
date
end_time=`date +%s`
runtime=$((end_time-start_time))
echo ${runtime}
