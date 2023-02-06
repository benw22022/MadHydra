export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "views LCG_102b_ATLAS_2 x86_64-centos7-gcc11-opt"
asetup 21.6.1, AthGeneration
source ${LCG_RELEASE_BASE}/LCG_88/MCGenerators/rivet/3.1.2/${LCG_PLATFORM}/rivetenv.sh
rivet-build -r /afs/cern.ch/work/b/bewilson/doubly_charged_higgs_generation/TypeII_joboptions/rivet/routines/pp_dppjj_llvvjj.cc |& tee compile.log
