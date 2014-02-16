#!/bin/bash
#################################
# PSI Tier-3 example batch Job  #
#################################

export DPM_HOST=se1.particles.ipm.ac.ir
export DPNS_HOST=se1.particles.ipm.ac.ir


##### CONFIGURATION ##############################################
# Output files to be copied back to the User Interface
# (the file path must be given relative to the working directory)
OUTFILES=

# Output files to be copied to the SE
# (as above the file path must be given relative to the working directory)
SEOUTFILES="output_$myVar.root"

SOURCEFILES=sourcefile

#
# By default, the files will be copied to $USER_SRM_HOME/$username/$JOBDIR,
# but here you can define the subdirectory under your SE storage path
# to which files will be copied (uncomment line)
#SEUSERSUBDIR="mytestsubdir/somedir"
#
# User's CMS hypernews name (needed for user's SE storage home path
# USER_SRM_HOME below)

HN_NAME=hpname

# set DBG=1 for additional debug output in the job report files
# DBG=2 will also give detailed output on SRM operations
DBG=0

#### The following configurations you should not need to change
# The SE's user home area (SRMv2 URL)
#USER_SRM_HOME="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/" saeid
USER_SRM_HOME="srm://se1.particles.ipm.ac.ir:8446/srm/managerv2?SFN=/dpm/particles.ipm.ac.ir/home/cms/store/user/"

# Top working directory on worker node's local disk. The batch
# job working directory will be created below this
#TOPWORKDIR=/scratch/`whoami`
TOPWORKDIR=/dataLOCAL/`whoami`

# Basename of job sandbox (job workdir will be $TOPWORKDIR/$JOBDIR)
JOBDIR=jobdir

SEUSERSUBDIR=srmdir

exe=exefile
##################################################################

ARGUMENTS=argumentsname

RELEASE=therelease

############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N example_job

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#  #$ -o /shome/username/mydir/
#  #$ -e /shome/username/mydir/

##################################################################



##### MONITORING/DEBUG INFORMATION ###############################
DATE_START=`date +%s`
echo "Job started at " `date`
cat <<EOF
################################################################
## QUEUEING SYSTEM SETTINGS:
HOME=$HOME
USER=$USER
JOB_ID=$JOB_ID
JOB_NAME=$JOB_NAME
HOSTNAME=$HOSTNAME
TASK_ID=$TASK_ID
QUEUE=$QUEUE

EOF

echo "################################################################"

if test 0"$DBG" -gt 0; then
   echo "######## Environment Variables ##########"
   env
   echo "################################################################"
fi


##### SET UP WORKDIR AND ENVIRONMENT ######################################
STARTDIR=`pwd`
WORKDIR=$TOPWORKDIR/$JOBDIR
RESULTDIR=$STARTDIR/$JOBDIR
if test x"$SEUSERSUBDIR" = x; then
   SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$JOBDIR
else
   SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$SEUSERSUBDIR
fi
if test -e "$WORKDIR"; then
   echo "WARNING: WORKDIR ($WORKDIR) already exists! Removing it..." >&2
   rm -rf $WOKDIR
fi
mkdir -p $WORKDIR
if test ! -d "$WORKDIR"; then
   echo "ERROR: Failed to create workdir ($WORKDIR)! Aborting..." >&2
   exit 1
fi

cd $WORKDIR
cat <<EOF
################################################################
## JOB SETTINGS:
STARTDIR=$STARTDIR
WORKDIR=$WORKDIR
RESULTDIR=$RESULTDIR
SERESULTDIR=$SERESULTDIR
EOF

###########################################################################
## YOUR FUNCTIONALITY CODE GOES HERE
# set up CMS environment

export SCRAM_ARCH=slc5_amd64_gcc462

source $VO_CMS_SW_DIR/cmsset_default.sh

cd $RELEASE

eval `scramv1 runtime -sh`

cd $WORKDIR


#saeid:: we do use dpm instead of dcap, needs to be checked
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
rootoutfiles=""
counter=0
for f in $SOURCEFILES; do
  output=${SEOUTFILES%.*}"_"$counter.root
  echo Running $exe $ARGUMENTS -o $output $f ...
  cp /dataLOCAL/MT2Tau/libshift.so.2.1  .
  export LD_LIBRARY_PATH=${pwd}:${LD_LIBRARY_PATH}
  $exe $ARGUMENTS -o $output $f 
  if test -e $output; then rootoutfiles=$rootoutfiles" "$output; fi
  let counter++
  CURRDATE=`date +%s`
  RUNTIME=$(echo "($CURRDATE - $DATE_START)/(60)" | bc -l)
  echo "-----> running since $RUNTIME min <-------------------------"
  RUNMIN=$(echo $RUNTIME | awk -F . '{print $myVar}')
  echo $RUNMIN
  if [ "$RUNMIN" -gt "55" ]; then echo "WARNING: runtime hits 1 hour.."; fi
done
if test x"$rootoutfiles" != x; then
  echo Running hadd $SEOUTFILES $rootoutfiles
  hadd $SEOUTFILES $rootoutfiles
fi

#### RETRIEVAL OF OUTPUT FILES AND CLEANING UP ############################
cd $WORKDIR
if test 0"$DBG" -gt 0; then
    echo "########################################################"
    echo "############# Working directory contents ###############"
    echo "pwd: " `pwd`
    ls -Rl
    echo "########################################################"
    echo "YOUR OUTPUT WILL BE MOVED TO $RESULTDIR"
    echo "########################################################"
fi

if test x"$OUTFILES" != x; then
   mkdir -p $RESULTDIR
   if test ! -e "$RESULTDIR"; then
          echo "ERROR: Failed to create $RESULTDIR ...Aborting..." >&2
          exit 1
   fi
   for n in $OUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          cp -a $WORKDIR/$n $RESULTDIR/$n
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORDFIR/$n to $RESULTDIR/$n" >&2
          fi
   fi
   done
fi

if test -e $SEOUTFILES; then
   if test 0"$DBG" -ge 2; then
      srmdebug="-debug"
   fi
   for n in $SEOUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          lcg-cp -b -D srmv2 file:///$WORKDIR/$n $SERESULTDIR/$n
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORKDIR/$n to $SERESULTDIR/$n" >&2
          fi
   fi
   done
fi

echo "Cleaning up $WORKDIR"
rm -rf $WORKDIR

###########################################################################
DATE_END=`date +%s`
RUNTIME=$((DATE_END-DATE_START))
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0

