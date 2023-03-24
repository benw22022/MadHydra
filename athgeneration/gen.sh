#!/bin/bash

DSID=$1
OUTPUT=$2
NEVENTS=$3

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
asetup AthGeneration, 21.6.90    # Need to use an older version of AthGeneration since JO was written in python2
export ATHENA_PROC_NUMBER=12

Gen_tf.py --ignorePatterns="attempt to add a duplicate" --ecmEnergy=13000. --maxEvents=$3 --firstEvent=1 --randomSeed=1234567 --outputEVNTFile=$2 --jobConfig=$1 --asetup=''