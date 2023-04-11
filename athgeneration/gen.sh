#!/bin/bash

DSID=$1 # "/home/bewilson/MadHydra/athgeneration/JO.py" #$1
OUTPUT=$2
NEVENTS=$3

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS="${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" 
setupATLAS -2
asetup AthGeneration, 21.6.40   # Need to use an older version of AthGeneration since JO was written in python2
export ATHENA_PROC_NUMBER=12

Gen_tf.py --ecmEnergy 13000 --maxEvents ${NEVENTS} --firstEvent 1 --randomSeed 1234556 --outputEVNTFile ${2}.EVNT.root} --jobConfig $1