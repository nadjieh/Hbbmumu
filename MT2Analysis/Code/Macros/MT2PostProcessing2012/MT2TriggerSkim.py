#!/usr/bin/python
##################################################################
# Script to copy MT2trees from T3 SE to the user's home          #
# -> MT2trees can be skimmed according to specified Triggers     #
#                                                                #
# Mario Masciovecchio                      February 7th, 2013    #
# from Pascal Nef's MT2PostProcessing.py                         #
# ################################################################
import sys, os, commands, shlex
from optparse import OptionParser
from sys import argv,exit

def MT2TriggerSkim():
	# parse arguments
	usage = """
	usage: %prog [options] file"
	
	This script will copy MT2trees from the T3 SE to the user's home directory 
	using the T3 batch system.
	-> use the --skim option to skim the MT2trees according to your 
	   preferred selection cuts. 

	Argument:
	file=lfn of the sample on the T3 SE
	-> in case the specified path contains subdirectories, the script will process
	   all the subdirectories contraining MT2trees (here you can wildcard subdirectories)

	Options:
	 --skim=shlib:    set path to shlib of MT2trees that was used to generate the MT2trees. 
	                  make sure the shlib corresponds to the correct MT2tag!
	                  in case this option is set, the MT2trees will be skimmed according 
		          to the cuts specified in 'MT2treeSkimming.C'
	 --prefix1=PREFIX: prefix for skimmed MT2trees, default is "TriggerSkimmed_MET"
	 --prefix2=PREFIX: prefix for skimmed MT2trees, default is "TriggerSkimmed_HTMHT"

	 --verbose:       set this option to get more prints
	 --script:        script name: default is MT2treeSkimming.C
	 --dryRun:        if this option is set, no jobs to the T3 batch system
	                  will be submitted. 

	Example:
	./MT2TriggerSkim.py --skim=~/MT2Analysis/Code/MT2AnalysisCode/MT2Code/shlib/libDiLeptonAnalysis.so --prefix=MyTriggerSkim [trigger1] [file1] [trigger2] [file2]
	with 
	[file1]=/store/user/mmasciov/SUSY/MassTrees/MT2_V02-02-00/MET/MET-Run2012D-PromptReco-v1-AOD/
	[file2]=/store/user/mmasciov/SUSY/MassTrees/MT2_V02-02-00/HTMHT/HTMHT-Run2012D-PromptReco-v1-AOD/
	[trigger1]='(trigger.HLT_PFMET150_v2||trigger.HLT_PFMET150_v3||trigger.HLT_PFMET150_v4||trigger.HLT_PFMET150_v5||trigger.HLT_PFMET150_v6||trigger.HLT_PFMET150_v7)'
	[trigger2]='((trigger.HLT_PFHT350_PFMET100_v4||trigger.HLT_PFHT350_PFMET100_v5||trigger.HLT_PFHT350_PFMET100_v6||trigger.HLT_PFHT350_PFMET100_v7||trigger.HLT_PFNoPUHT350_PFMET100_v1||trigger.HLT_PFNoPUHT350_PFMET100_v3||trigger.HLT_PFNoPUHT350_PFMET100_v4)&&!(trigger.HLT_PFMET150_v2||trigger.HLT_PFMET150_v3||trigger.HLT_PFMET150_v4||trigger.HLT_PFMET150_v5||trigger.HLT_PFMET150_v6||trigger.HLT_PFMET150_v7))'
	"""

	parser = OptionParser(usage)

	parser.add_option("--skim"   ,dest="shlib",      help="shlib: path to shlib used to generate the MT2trees")
	parser.add_option("--MT2tag" ,dest="tag",        help="MT2tag: to set the tag manually")
	parser.add_option("--prefix1" ,dest="prefix1",     help="prefix1: default is 'TriggerSkimmed_MET'")
	parser.add_option("--prefix2" ,dest="prefix2",     help="prefix2: default is 'TriggerSkimmed_HTMHT'")	
	parser.add_option("--script" ,dest="script",     help="script: default is 'MT2treeTriggerSkimming.C'")
	parser.add_option("--verbose",dest="verbose",    help="set verbose", action="store_true")
	parser.add_option("--dryRun", dest="dryRun",     help="dry Run: do not do anything", action="store_true")
	#parser.add_option("--rootFile", dest="rootFile", help="rootFile: choose root files to process (accepts wildcards)")
	#parser.add_option("--trigger1" ,dest="trigger1",     help="trigger1: default is 'HLT_PFMET150'")
	#parser.add_option("--trigger2" ,dest="trigger2",     help="trigger2: default is 'HLT_PFHT350_PFMET100'")

	global options, args, MT2tag
	(options, args) = parser.parse_args()

	if(options.verbose): print options
	
	if(len(args)!=4):
	        print "exactly four arguments must be given!"
		print "try --help option"
		exit(-1)

	trigger1 = args[0]
	FILE1   = args[1]

	trigger2 = args[2]
	FILE2   = args[3]

	if not os.path.exists("logs"):
		print "creating dir 'logs' for log-files"
		os.mkdir("logs")

	if not options.tag==None:
		MT2tag = options.tag
	else:
		MT2tag = FILE1[FILE1.find("/MT2_V")+1:FILE1.find("/",FILE1.find("/MT2_V")+1)]
		if(len(MT2tag)!=13):
			print "error parsing MT2_tag from argument"
			exit(-1)
		else:
			print "setting MT2tag= " +MT2tag
	if options.script==None and not options.shlib==None:
		print "using script MT2treeTriggerSkimming.C"
	
	#FILE 1: MET dataset
	WILDCARD = ""
	if (FILE1.find("*")!=-1):
	        indx = FILE1.rfind("/",0,FILE1.find("*"))
		WILDCARD = FILE1[indx+1:].replace("*",".*")
		FIL1E = FILE1[0:indx]
	t3_se_dir       ="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/"
	if (WILDCARD!=""):
	        command   ="srmls "+ t3_se_dir + FILE1+ " | grep -v .root | grep '" + WILDCARD + "' | awk '{print $2}' | awk -F trivcat '{print $2}'"
	else:
	        command   ="srmls "+ t3_se_dir + FILE1+ " | grep -v .root | awk '{print $2}' | awk -F trivcat '{print $2}'"
	status, output  = commands.getstatusoutput(command)
	if options.verbose:
	        print command
		print output  
	subdirs =output.splitlines()
	if(len(subdirs) >1): subdirs.pop(0)  # remove base dir from list in case there are sub-directories 
	subdirs.sort()

	for subdir in subdirs:
	        #if (options.rootFile!=None):
		        #command   ="srmls "+ t3_se_dir + subdir +" | grep .root | grep '" + options.rootFile.replace('*','.*') + "' | awk '{print $2}' | awk -F trivcat '{print $2}' | sort "
		#else:
		command   ="srmls "+ t3_se_dir + subdir +" | grep .root | awk '{print $2}' | awk -F trivcat '{print $2}' | sort "
		status, files  = commands.getstatusoutput(command)
		if(options.verbose):
			print files
			print command
			print ""
		filelist =""
		for file in files.splitlines():
			filelist += file + " " 
		os.system("rm -rf "+logfile(subdir))
		##cmd = "qsub -q short.q -N TriggerSkim_"+sample(subdir)+" -o "+logfile(subdir)+" -j y  ./MT2TriggerSkim.csh "
		
		cmd = "qsub -q all.q -N TriggerSkim_"+sample(subdir)+" -o "+logfile(subdir)+"  -j y ./MT2TriggerSkim.csh "
		cmd = cmd + " "+ MT2tag + " "+ sample(subdir)

		if(options.shlib!= None):
			cmd = cmd + " "+ options.shlib
		else:
			cmd =" "
		if(options.prefix1!=None):
			cmd = cmd+ " " + options.prefix1
		else:
			cmd = cmd+ " " + "TriggerSkimmed_MET"
		if(options.script!=None):
			cmd = cmd+ " " + options.script
		else:
			cmd = cmd+ " " + "MT2treeTriggerSkimming.C"

		cmd = cmd +" '"+trigger1+"' "+ filelist
		print cmd
		
		if not options.dryRun:
			os.system(cmd)
		print "---------------------------------------------------------------------"

	#FILE 2: HTMHT
	if not options.tag==None:
		MT2tag = options.tag
	else:
		MT2tag = FILE2[FILE2.find("/MT2_V")+1:FILE2.find("/",FILE2.find("/MT2_V")+1)]
		if(len(MT2tag)!=13):
			print "error parsing MT2_tag from argument"
			exit(-1)
		else:
			print "setting MT2tag= " +MT2tag
	if options.script==None and not options.shlib==None:
		print "using script MT2treeTriggerSkimming.C"

	if (FILE2.find("*")!=-1):
	        indx = FILE2.rfind("/",0,FILE2.find("*"))
		WILDCARD = FILE2[indx+1:].replace("*",".*")
		FIL1E = FILE2[0:indx]
	t3_se_dir       ="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/"
	if (WILDCARD!=""):
	        command   ="srmls "+ t3_se_dir + FILE2+ " | grep -v .root | grep '" + WILDCARD + "' | awk '{print $2}' | awk -F trivcat '{print $2}'"
	else:
	        command   ="srmls "+ t3_se_dir + FILE2+ " | grep -v .root | awk '{print $2}' | awk -F trivcat '{print $2}'"
	status, output  = commands.getstatusoutput(command)
	if options.verbose:
	        print command
		print output  
	subdirs =output.splitlines()
	if(len(subdirs) >1): subdirs.pop(0)  # remove base dir from list in case there are sub-directories 
	subdirs.sort()

	for subdir in subdirs:
		command   ="srmls "+ t3_se_dir + subdir +" | grep .root | awk '{print $2}' | awk -F trivcat '{print $2}' | sort "
		status, files  = commands.getstatusoutput(command)
		if(options.verbose):
			print files
			print command
			print ""
		filelist =""
		for file in files.splitlines():
			filelist += file + " " 
		os.system("rm -rf "+logfile(subdir))
		#cmd = "qsub -q short.q -N TriggerSkim_"+sample(subdir)+" -o "+logfile(subdir)+" -j y  ./MT2TriggerSkim.csh "
		cmd = "qsub -q all.q -N TriggerSkim_" + sample(subdir) + " -o "+logfile(subdir)+" ./MT2TriggerSkim.csh "
		cmd = cmd + " "+ MT2tag + " "+ sample(subdir)

		if(options.shlib!= None):
			cmd = cmd+ " "+ options.shlib
		else:
			cmd =" "
		if(options.prefix2!=None):
			cmd = cmd+ " " + options.prefix2
		else:
			cmd = cmd+ " " + "TriggerSkimmed_HTMHT"
		if(options.script!=None):
			cmd = cmd+ " " + options.script
		else:
			cmd = cmd+ " " + "MT2treeTriggerSkimming.C"
			
		cmd = cmd + " '"+trigger2+"' "+filelist
		print cmd
		
		if not options.dryRun:
			os.system(cmd)
		print "---------------------------------------------------------------------"

def logfile(subdir):
	logname = subdir.split(MT2tag)[1]
	logname = logname.replace("/","_")
	logname = logname.strip("_")
	#logname = "logs/"+logname+".log"

        logname = "logs/"+logname+"_TriggerSkim.log"

	print "   THIS IS THE LOGNAME: " + logname
	return logname

def sample(subdir):
	sample = subdir.split(MT2tag)[1]
	sample = sample.strip("/")
	sample = sample.split("/")[1]
	sample = sample.strip("/")
	print "   THIS IS THE SAMPLE: " + sample
	print " This is Subdir: " + subdir
	return sample

if __name__ == "__main__":
	MT2TriggerSkim()
