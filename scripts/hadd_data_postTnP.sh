#!/usr/bin/sh
#import os

for pathname_full in /uscms/home/amalbert/nobackup/CMSSW_10_6_20/src/run3_llp_analyzer/Merged_Cache_InputLists/*; do
	pathname=$(basename $pathname_full)
	echo $pathname
	# if [[ $pathname == *"2022"* ]]; then
    #     	filebase=$(echo $pathname | cut -d. -f1)
	# 	hadd ${filebase}.root `xrdfsls -u /store/group/lpclonglived/amalbert/Data_MC_Comp_TnP/results_from_cache_noSkim/Data_try2/2022/$filebase | grep '\.root'` 
    #     	xrdcp ${filebase}.root root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/Data_MC_Comp_TnP/results_from_cache_noSkim/Data_try2/2022_Merged
	# 	rm ${filebase}.root
    # 	fi

    if [[ $pathname == *"2023C"* ]]; then
        filebase=$(echo $pathname | cut -d. -f1)
        hadd ${filebase}.root `xrdfsls -u /store/group/lpclonglived/amalbert/Data_MC_Comp_TnP/results_from_cache_noSkim/Data_noClusters/2023/$filebase | grep '\.root'` 
        xrdcp ${filebase}.root root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/Data_MC_Comp_TnP/results_from_cache_noSkim/Data_noClusters/2023_Merged
	rm ${filebase}.root
    fi
done

