#!/bin/bash

DSID=$1 # "/home/bewilson/MadHydra/athgeneration/JO.py" #$1
OUTPUT=$2
NEVENTS=$3

echo \n \n \n $1 \n \n \n

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh -2
asetup AthGeneration,  21.6.57   # Need to use an older version of AthGeneration since JO was written in python2
export ATHENA_PROC_NUMBER=12

cat ../../athgeneration/JOB/HeavyN_emu_1TeV.py

Gen_tf.py --ignorePatterns="attempt to add a duplicate" --ecmEnergy=13000 --maxEvents=${NEVENTS} --firstEvent=1 --randomSeed=123456 --outputEVNTFile=${2}.EVNT.root --jobConfig="../../athgeneration/JOB/HeavyN_emu_1TeV.py" --asetup=''