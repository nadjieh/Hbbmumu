#!/usr/bin/python

import sys, os, commands

serious = ['"core dump"',
	   '"segmentation violation"', 
	   '"Killed"', 
	   '"probably not closed"', 
	   '"credential remaining lifetime"', 
	   '"SRM_FAILURE"',
	   '"ERROR: there is something wrong"',
	   '"ERROR: Failed to copy"'] 

timing = ['"WARNING: runtime hits 1 hour.."']

to_check = ['"skip event"', 
	    '"WARNING: crazy HCAL event"', 
	    '"ERROR: cannot open dead cell file"', 
	    '"has NAN-Pt Objects!! Skipping Event"', 
	    '"WARNING: MT2 calculation failed"', 
	    '"WARNING: MT2 minimizeDHT calculation failed"', 
#	    '"WARNING: Event with Jets with negative JEC"'
	    ]


os.system('echo "************SERIOUS PROBLEMS*******************"')
for i in serious:
	cmd = "grep -ri --exclude='myScript*.sh' " + i + " ."
	print cmd
	os.system(cmd)

os.system('echo "************ TIMING ***************************"')
for i in timing:
	cmd = "grep -ri --exclude='myScript*.sh' " + i + " ."
	print cmd
	os.system(cmd)

os.system('echo "************WARNINGS***************************"')
for i in to_check:
	cmd = "grep -ri " + i + " ."
	print cmd
	os.system(cmd)
