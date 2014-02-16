#!/usr/bin/python

### WATCH OUT: NEEDS TO BE LAUNCHED AS python runManager.py configuration.cfg and NOT ./runManager configuration.cfg

# System libraries to interact with the shell
import sys
import os
from os import popen
import commands
import time

from ROOT import TFile, TTree, gSystem, gROOT

from DBSAPI.dbsApi import DbsApi
#from DBSAPI.dbsException import *
#from DBSAPI.dbsApiException import *


def root_logon():
  gSystem.SetIncludePath(" -I./include/ -I../TESCO/include/");
  #gSystem->Load("/public/V_CMSSW64/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/slc5_amd64_gcc462/libFWCoreFWLite.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreFWLite.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreUtilities.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsCommon.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsFWLite.so");

  gSystem.Load("libPhysics");
  gSystem.Load("shlib/libDiLeptonAnalysis.so");
  gROOT.ProcessLine(".x SetStyle_PRD.C");

def showMessage(message):
  print "[runManager] " + message


def process(dataset,how_to_list='file'): #command,python,file
  dpmPath = "rfio:///dpm/particles.ipm.ac.ir/home/cms/"
  dbsUrl = 'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet'
  list = [] # Initialize file list

  # Check if cache file exists and if we should use it
  cacheFile = '.'+dataset.replace('/','_').replace('?','_')
  allmyfiles=[]
  if how_to_list == 'file' and os.path.isfile(cacheFile):
      print "Reading from cached file"
      f = open(cacheFile)
      showMessage('Reading files pertaining to '+dataset+'from cache file '+cacheFile)
      allmyfiles = f.readlines()
      f.close()
  elif how_to_list == 'python':
      print "Trying Good PATH"
      api = DbsApi({'url':dbsUrl})
      for afile in api.listDatasetFiles(datasetPath=dataset):
        allmyfiles.append( afile["LogicalFileName"] + " " )
        print afile["LogicalFileName"]
      f = open(cacheFile,'w')
      f.write(''.join(allmyfiles))
  else:
      # If not: rebuild the list
      showMessage("Going to fetch list of files pertaining to "+str(dataset))
      command = 'dbs search --query="find file where dataset='+dataset+'" --noheader --url='+dbsUrl
      print command
      theList = os.popen(command)
      list = theList.readlines()
      for li in list:
              allmyfiles.append(li)
      f = open(cacheFile,'w')
      f.write(''.join(allmyfiles))
  
  numberOfFiles = len(allmyfiles)

  correctList = [];
  for fileName in allmyfiles:
    auxiliar = dpmPath + fileName
    auxiliar = auxiliar[0:len(auxiliar)-1]
    correctList.append(auxiliar)



  root_logon()

  
  
  N = 0
  NEntriesTotal = 0
  f2 = open(cacheFile+"_summary.csv",'w')
  for fileName in correctList:
    file = TFile.Open(fileName)
    t = file.Get("Events")
    N = N + 1
    NEntriesTotal = NEntriesTotal + t.GetEntries()
    print str(N) + " " + fileName + ", NEntries:" + str(t.GetEntries()) + ", NTotal:" + str(NEntriesTotal)
    print >> f2 ,"%(N)d,%(filename)s,%(NEntries)d,%(NTotal)d" % {"N":N , "filename":fileName , "NEntries":t.GetEntries(), "NTotal":NEntriesTotal}
#    f2.flush()
    file.Close()
    
#process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012A-13Jul2012-v1-61aef9e3bbf721ce333708258d830dbb/USER")          DONE
#process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012A-recover-06Aug2012-v1-61aef9e3bbf721ce333708258d830dbb/USER","python")DONE  
#process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012B-13Jul2012-v1-61aef9e3bbf721ce333708258d830dbb/USER","python") DONE
#process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012C-24Aug2012-v1-1c5ff4400f073c7d50e40dc88ea10d23/USER","python") DONE
#process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012C-PromptReco-v2-2c909eaf42bba0fead9d3c5a78eda5f3/USER") DONE
#process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012D-PromptReco-v1-2c909eaf42bba0fead9d3c5a78eda5f3/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012B-05Nov2012-v2-61aef9e3bbf721ce333708258d830dbb/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012C-part1-05Nov2012-v2-1243ed21a886eafc37154ce6d942024f/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012C-part2-05Nov2012-v2-FirstHalf-1243ed21a886eafc37154ce6d942024f/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012C-part2-05Nov2012-v2-SecondHalf-1243ed21a886eafc37154ce6d942024f/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012D-part1-10Dec2012-v1-FirstHalf-913e379464056a6b529e5016356a02b3/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012D-part1-10Dec2012-v1-SecondHalf-913e379464056a6b529e5016356a02b3/USER") DONE


#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012B-05Nov2012-v2-RemainingRereco-61aef9e3bbf721ce333708258d830dbb/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012C-part2-05Nov2012-v2-RemainingRereco-1243ed21a886eafc37154ce6d942024f/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012D-part1-10Dec2012-v1-RemainingRereco-913e379464056a6b529e5016356a02b3/USER") DONE
#process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012D-part2-17Jan2013-v1-RemainingRereco-b25f9a28f42fdaa751ccabd67a60075a/USER") DONE
process("/MultiJet1Parked/paktinat-V03-09-01_MultiJet1Parked-Run2012D-part2-PixelRecover-17Jan2013-v1-b25f9a28f42fdaa751ccabd67a60075a/USER") #DONE
