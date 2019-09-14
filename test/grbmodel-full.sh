#!/bin/sh

##########################################################
# Shell script to launch a toy model GRB analysis job    #
# to the batch system of the IIHE IceCube cluster.       #
#                                                        #
# The actual C++ code is contained in the file           #
# grbmodel.cxx                                           #
#                                                        #
# To launch the job(s) issue the command ($ = prompt)    #
#                                                        #
# $ source grbmodel.sh                                   #
#                                                        #
# The GCN data files are specified in this script.       #
# The log file will have the job id name with the        #
# extension ".log" and the produced output file will     #
# obtain the extension ".root".                          #
#                                                        #
# In the C++ code itself just some generic filenames     #
# will be used which will only exist as temporary files  #
# on the local disks of the cluster during processing.   #
#                                                        #
# All output and log files will be returned into the     #
# (working) directory from which the job was launched.   #
#                                                        #
# Nick van Eijndhoven, 22-sep-2014, IIHE-VUB, Brussels   #
##########################################################

### The (generic) GRB identifier
export GRB=GRB

### The Queue selection for the batch system (localgrid (=any), express, highmem or gpu)
export QUEUE=express

###
### Create and submit a shell script "process-job.sh" to run the programme
###

### Construction of the name of the job and log and output file basename for this job
export BASENAME=grbmodel

###
### Creation of the corresponding shell script for this subrun
###

cat >process-job.sh <<*EOR
#################################################################################
### This is an automatically generated shell script for PBS batch processing. ###
### The job may be submitted via the following command : ($ = prompt)         ###
### $ qsub process-job.sh                                                     ###
#################################################################################
#!/bin/sh

### Combine the error and output files for this job
#PBS -j oe

### Jobname specification
#PBS -N $BASENAME

### The log and output filename base for this subrun
export FILENAME=$BASENAME-\$PBS_JOBID

### The working directory at the host from which the job was submitted
export HOSTDIR=\$PBS_O_WORKDIR

### The name of the batch slave
export SLAVE=\$HOSTNAME

### The working directory at the batch slave
export SLAVEDIR=/scratch/\$PBS_JOBID

### Go to the batch machine working directory and clear it
cd \$SLAVEDIR
rm -rf *

### Import all necessary files
cp \$HOSTDIR/NcGRB-full.cxx .
cp \$HOSTDIR/grbmodel-full.cc .
cp /data/IC86/exp/GRB/GCN/GRB-t90-z-Swift.txt .
cp /data/IC86/exp/GRB/GCN/GRB-t90-Fermi.txt .
cp /data/IC86/exp/GRB/GCN/GRB-IC86*.root .

echo Contents of the working directory \$SLAVEDIR on \$SLAVE for job \$PBS_JOBNAME
ls -lt

### Define the log output file for this job
LOGFILE=\$HOSTDIR/\$FILENAME.log
 
### Initialise the ROOT/NCFS environment
source /ice3/software/iihe/ncfs.sh

### Execute the ROOT macro
root -b -q grbmodel-full.cc >\$LOGFILE

### Copy the produced output file to the corresponding host directory
cp output.root \$HOSTDIR/\$FILENAME.root

### Clear the batch machine working directory
cd \$SLAVEDIR
rm -rf *

### End of the PBS job processing script(s) ### 
*EOR

### Submit the shell script for this subrun processing
qsub -q $QUEUE process-job.sh

### End of job(s) submission ###
 