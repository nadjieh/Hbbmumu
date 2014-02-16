#!/usr/bin/python
##################################################################
# Script to copy MT2trees from T3 SE to the user's home          #
# -> MT2trees can be skimmed according to specified cuts.        #
#                                                                #  
# Pascal Nef                             September 18th, 2011    #  
# ################################################################
import sys, os, commands, shlex
from optparse import OptionParser
from sys import argv,exit

def MT2PostProcessing():
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
	 --prefix=PREFIX: prefix for skimmed MT2trees, default is "skimmed" 
	 --verbose:       set this option to get more prints
	 --script:        script name: default is MT2treeSkimming.C
	 --dryRun:        if this option is set, no jobs to the T3 batch system
	                  will be submitted. 
	 --rootFile:      process only given root files (accepts wildcards)

	Example:
	./MT2PostProcessing.py --skim=~/MT2Analysis/Code/MT2AnalysisCode_V01-00-00/MT2Code/shlib/libDiLeptonAnalysis.so --prefix=skimmed_highMT2 [file]  --rootFile=myfile*.root
	with [file]=/store/user//pnef/SUSY/MassTrees/MT2_V01-00-00/20110915_MSSM_MC_nocuts/QCD*/

	"""
	parser = OptionParser(usage)

	parser.add_option("--skim"   ,dest="shlib",      help="shlib: path to shlib used to generate the MT2trees")
	parser.add_option("--MT2tag" ,dest="tag",        help="MT2tag: to set the tag manually")
	parser.add_option("--prefix" ,dest="prefix",     help="prefix: default is 'skimmed'")
	parser.add_option("--script" ,dest="script",     help="script: default is 'MT2treeSkimming.C'")
	parser.add_option("--verbose",dest="verbose",    help="set verbose", action="store_true")
	parser.add_option("--dryRun", dest="dryRun",     help="dry Run: do not do anything", action="store_true")
	parser.add_option("--rootFile", dest="rootFile", help="rootFile: choose root files to process (accepts wildcards)")

	global options, args, MT2tag
	(options, args) = parser.parse_args()

	if(options.verbose): print options
	
	if(len(args)!=1):
	        print "exactly one argument must be given!"
		print "try --help option"
		exit(-1)

	FILE   = args[0]

	if not os.path.exists("logs"):
		print "creating dir 'logs' for log-files"
		os.mkdir("logs")

	if not options.tag==None:
		MT2tag = options.tag
	else:
		MT2tag = FILE[FILE.find("/MT2_V")+1:FILE.find("/",FILE.find("/MT2_V")+1)]
		if(len(MT2tag)!=13):
			print "error parsing MT2_tag from argument"
			exit(-1)
		else:
			print "setting MT2tag= " +MT2tag
	if options.script==None and not options.shlib==None:
		print "using script MT2treeSkimming.C"
	
	# get list of files
	WILDCARD = ""
	if (FILE.find("*")!=-1):
	        indx = FILE.rfind("/",0,FILE.find("*"))
		WILDCARD = FILE[indx+1:].replace("*",".*")
		FILE = FILE[0:indx]
	t3_se_dir       ="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/"
	if (WILDCARD!=""):
	        command   ="srmls "+ t3_se_dir + FILE+ " | grep -v .root | grep '" + WILDCARD + "' | awk '{print $2}' | awk -F trivcat '{print $2}'"
	else:
	        command   ="srmls "+ t3_se_dir + FILE+ " | grep -v .root | awk '{print $2}' | awk -F trivcat '{print $2}'"
	status, output  = commands.getstatusoutput(command)
	if options.verbose:
	        print command
		print output  
	subdirs =output.splitlines()
	if(len(subdirs) >1): subdirs.pop(0)  # remove base dir from list in case there are sub-directories 
	subdirs.sort()

	for subdir in subdirs:
	        if (options.rootFile!=None):
		        command   ="srmls "+ t3_se_dir + subdir +" | grep .root | grep '" + options.rootFile.replace('*','.*') + "' | awk '{print $2}' | awk -F trivcat '{print $2}' | sort "
		else:
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
		#cmd = "qsub -q short.q -N MT2_"+sample(subdir)+" -o "+logfile(subdir)+" -j y  ./MT2PostProcessing.csh "
		#MM
		cmd = "qsub -q all.q -N MT2_"+sample(subdir)+" -o "+logfile(subdir)+" -j y  ./MT2PostProcessing.csh "
		#MM
		cmd = cmd + " "+ MT2tag + " "  + production(subdir) +  " "  + sample(subdir) 
		if(options.shlib!= None):
			cmd = cmd+ " " + options.shlib
		else:
			cmd = cmd+ " " + "none"
		if(options.prefix!=None):
			cmd = cmd+ " " + options.prefix
		else:
			cmd = cmd+ " " + "skimmed"
		if(options.script!=None):
			cmd = cmd+ " " + options.script
		else:
			cmd = cmd+ " " + "MT2treeSkimming.C"
		cmd = cmd + " "+ filelist
		print cmd
		if not options.dryRun:
			os.system(cmd)
		print "---------------------------------------------------------------------"


def logfile(subdir):
	logname = subdir.split(MT2tag)[1]
	logname = logname.replace("/","_")
	logname = logname.strip("_")
	logname = "logs/"+logname+".log"
##MM
        #logname = "logs/"+logname+"_HTcut.log"
##MM
	print "   THIS IS THE LOGNAME: " + logname
	return logname

def production(subdir):
	production = subdir.split(MT2tag)[1]
	production = production.strip("/")
	production = production.split("/")[0]
	print "   THIS IS THE PRODUCTION: " + production
	return production

def sample(subdir):
	#sample = subdir.split(production(subdir))[1]
	sample = subdir.split(MT2tag)[1]
	sample = sample.strip("/")
	sample = sample.split("/")[1]
	sample = sample.strip("/")
	print "   THIS IS THE SAMPLE: " + sample
	print " This is Subdir: " + subdir
	return sample

if __name__ == "__main__":
	MT2PostProcessing()
