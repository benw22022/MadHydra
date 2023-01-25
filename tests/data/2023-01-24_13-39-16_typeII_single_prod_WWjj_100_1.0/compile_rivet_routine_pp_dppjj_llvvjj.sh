echo hello
setupATLAS
lsetup root
asetup 21.6.1, AthGeneration
source ${LCG_RELEASE_BASE}/LCG_88/MCGenerators/rivet/3.1.2/${LCG_PLATFORM}/rivetenv.sh
rivet-build /afs/cern.ch/work/b/bewilson/doubly_charged_higgs_generation/TypeII_joboptions/rivet/routines/pp_dppjj_llvvjj.cc
