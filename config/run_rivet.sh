export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
asetup 23.6.1, AthGeneration
lsetup "views LCG_102b_ATLAS_2 x86_64-centos7-gcc11-opt"
rivet-build -r ROUTINE.cc |& tee compile.log
rivet --analysis=ROUTINE --pwd INPUTFILE |& tee run.log
