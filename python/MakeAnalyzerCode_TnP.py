#!/usr/bin/env python

import sys
#print("IN FILE")
if len(sys.argv) != 2:
    print "Please specify an analyzer name!"
    sys.exit()

analyzer = sys.argv[1]

inNames = ['include/AnalyzerTemplate_TnP.txt','src/RunAnalyzerTemplate.txt']
outNames = ['analyzers/'+analyzer+'.h','src/Run'+analyzer+'.cc']
#print("outNames:", outNames)
if analyzer.find("Run1") > 0:
	inNames = ['include/AnalyzerTemplateRun1.txt','src/RunAnalyzerTemplateRun1.txt']

if analyzer.find("UpgradeTiming") > 0:
	inNames = ['include/AnalyzerTemplateUpgradeTiming.txt','src/RunAnalyzerTemplateUpgradeTiming.txt']

for i in range(len(inNames)):
    with open(inNames[i]) as inF:
        with open(outNames[i],'w') as outF:
            for line in inF:
                outF.write( line.replace('%ANALYZER%',analyzer) )       
