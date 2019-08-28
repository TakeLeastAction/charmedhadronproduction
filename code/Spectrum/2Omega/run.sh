#!/bin/bash

# Receive Arguments from command lines
path=$1
workpath=$2 
user=$3
jobNode=`hostname`

# Create an uniquename for the run folder
unique_tag1=`basename $path | cut -d'.' -f1-2`
unique_tag2=`date +%s`
localdir="/tmp/${user}/mytest/${workpath}-${unique_tag1}-${unique_tag2}"
workdir=$path"/"$workpath

echo "..................................................................................."
echo "  Runing Jobs using host $jobNode in directory $workdir..."
echo "..................................................................................."

#### create local working directory and copy input files to local disk

mkdir -vp $localdir

cp -rf $path/$workpath/*    $localdir/

#### change to local directory and run jobs
cd $localdir

atlasSetup
localSetupGcc
localSetupROOT

### run job
#make
./col_xn > output.txt

#### copy output files back to current working directory
cp  *       $workdir

#### delete information in local disk
cd /tmp/$user/mytest
rm -r -f $localdir

########################################################

ExitStatus=$?
if [ $ExitStatus -ne 0 ]; then
  echo "e-id job failed with status $ExitStatus ..."
  exit $ExitStatus
fi

exit 0

