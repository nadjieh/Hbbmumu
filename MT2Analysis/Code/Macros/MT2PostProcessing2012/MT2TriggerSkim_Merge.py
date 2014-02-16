#!/usr/bin/python
##################################################################
# Script to copy MT2trees from T3 SE to the user's home          #
# -> MT2trees can be skimmed according to specified Triggers     #
#                                                                #
# Mario Masciovecchio                      February 8th, 2013    #
# ################################################################
import sys, os, commands, shlex
from optparse import OptionParser
from sys import argv,exit

def MT2TriggerSkim_Merge():
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

	 --prefix[1,2]=PREFIX: prefix for skimmed MT2trees, default is "TriggerSkimmed_MET" or "TriggerSkimmed_HTMHT"]
	 --verbose:       set this option to get more prints
	 --dryRun:        if this option is set, no jobs to the T3 batch system
	                  will be submitted. 

	Example:
	./MT2TriggerSkim_Merge.py [file1] [file2]
	with [file1]=/shome/mmasciov/MT2Analysis/MT2trees/MT2_V02-02-00/TriggerSkimmed_MET/MET-Run2012D-PromptReco-v1-AOD.root
	[file2]=/shome/mmasciov/MT2Analysis/MT2trees/MT2_V02-02-00/TriggerSkimmed_HTMHT/HTMHT-Run2012D-PromptReco-v1-AOD.root
	"""
	
	parser = OptionParser(usage)

	parser.add_option("--MT2tag1" ,dest="tag1",        help="MT2tag1: to set the tag manually")
	parser.add_option("--MT2tag2" ,dest="tag2",        help="MT2tag2: to set the tag manually")
	parser.add_option("--prefix1" ,dest="prefix1",     help="prefix1: default is 'TriggerSkimmed_MET'")
	parser.add_option("--prefix2" ,dest="prefix2",     help="prefix2: default is 'TriggerSkimmed_HTMHT'")	
	parser.add_option("--verbose",dest="verbose",    help="set verbose", action="store_true")
	parser.add_option("--dryRun", dest="dryRun",     help="dry Run: do not do anything", action="store_true")

	global options, args, MT2tag1, MT2tag2
	(options, args) = parser.parse_args()

	if(options.verbose): print options
	
	if(len(args)!=2):
	        print "exactly two arguments must be given!"
		print "try --help option"
		exit(-1)

	FILE1   = args[0]
	FILE2   = args[1]

	if not os.path.exists("logs"):
		print "creating dir 'logs' for log-files"
		os.mkdir("logs")

	#FILE 1: MET dataset
	if not options.tag1==None:
		MT2tag1 = options.tag
	else:
		MT2tag1 = FILE1[FILE1.find("/MT2_V")+1:FILE1.find("/",FILE1.find("/MT2_V")+1)]
		if(len(MT2tag1)!=13):
			print "error parsing MT2_tag1 from argument"
			exit(-1)
		else:
			print "setting MT2tag1= " +MT2tag1
	
	if(options.prefix1!=None):
		prefix1=options.prefix1
	else:
		prefix1="TriggerSkimmed_MET"

		
	Sample1= FILE1[FILE1.find(prefix1+"/")+len(prefix1+"/"):FILE1.find(".root",FILE1.find(prefix1+"/")+len(prefix1+"/"))]

	#FILE 2: HTMHT
	if not options.tag2==None:
		MT2tag2 = options.tag
	else:
		MT2tag2 = FILE2[FILE2.find("/MT2_V")+1:FILE2.find("/",FILE2.find("/MT2_V")+1)]
		if(len(MT2tag2)!=13):
			print "error parsing MT2_tag2 from argument"
			exit(-1)
		else:
			print "setting MT2tag2= " +MT2tag2

	if(options.prefix2!=None):
		prefix2=options.prefix2
	else:
		prefix2="TriggerSkimmed_HTMHT"

	Sample2= FILE2[FILE2.find(prefix2+"/")+len(prefix2+"/"):FILE2.find(".root",FILE2.find(prefix2+"/")+len(prefix2+"/"))]

	os.system("rm -rf "+logfile())

	#Merging:
	cmd_merge="qsub -q all.q -N TriggerSkim_Merge_" + Sample1+"_"+Sample2+ " -o "+logfile() + " ./MT2TriggerSkim_Merge.csh "
	cmd_merge= cmd_merge+ " " + MT2tag1+ " " +prefix1+ " " + Sample1+ " "+MT2tag2+ " "+prefix2+ " " +Sample2

	print cmd_merge

	print "Merging Skimmed Ntuples..."
	os.system(cmd_merge)
	print "---------------------------------------------------------------------"
	#End of Merging
	
def logfile():
	logname="logs/MergingAttempt.log"
	print "   THIS IS THE LOGNAME: " + logname
	return logname

if __name__ == "__main__":
	MT2TriggerSkim_Merge()
