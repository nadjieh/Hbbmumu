#!/bin/tcsh -f
##################################################################
# Script to skim MT2trees on T3 batch system                     #
# Mario Masciovecchio                       February 8h, 2013    # 
# ################################################################

# Avoid wild cards
set noglob
############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

##### MONITORING/DEBUG INFORMATION ###############################
set DATE_START=`date +%s`
echo "Job started at " `date`

echo "Very beginning... Here we are!"
### Check n agruments ###########################################
if ($#argv < 6) then
	echo "usage "$0" MT2tag production sample shlib prefix script file(s)"
	exit(-1)
endif

echo "After checking N arguments..."

set NARGS      = $#argv
set MT2TAGa    = $1
set PREFIXa    = $2
set SAMPLEa    = $3
set MT2TAGb    = $4
set PREFIXb    = $5
set SAMPLEb    = $6

set MERGEDsample = $3_$6.root
echo $MERGEDsample

set HOME=`pwd`
set OUTDIR=/shome/`whoami`/MT2Analysis/MT2trees/TriggerSkimmed/Merged/
set TOPWORKDIR=/scratch/`whoami`
set WORKDIR=/scratch/`whoami`/Merging/
mkdir -pv $WORKDIR


########## Set the CMSSW environment (for ROOT, mainly...)
source $VO_CMS_SW_DIR/cmsset_default.csh
setenv SCRAM_ARCH slc5_amd64_gcc462
cd /shome/mmasciov/CMSSW_5_3_6/src
cmsenv
# Fix problem with DCache and ROOT
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
cd -

echo "nargs " $NARGS

cd $WORKDIR
cp /shome/`whoami`/MT2Analysis/MT2trees/$1/$2/$3.root .
cp /shome/`whoami`/MT2Analysis/MT2trees/$4/$5/$6.root .

hadd $MERGEDsample $3.root $6.root

echo "copying merged skimmed MT2tree"
mkdir -pv $OUTDIR/
cp -v $MERGEDsample $OUTDIR/
	
# clean up ################################################################
rm -r $WORKDIR 
echo -------------------- ALL DONE --------------------------------------
###########################################################################
set DATE_END=`date +%s`
@ RUNTIME = $DATE_END - $DATE_START
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0

