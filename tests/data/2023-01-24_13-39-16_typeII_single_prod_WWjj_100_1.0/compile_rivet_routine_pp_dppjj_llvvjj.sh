export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
asetup 23.6.1, AthGeneration
lsetup "views LCG_102b_ATLAS_2 x86_64-centos7-gcc11-opt"
rivet-build -r /afs/cern.ch/work/b/bewilson/doubly_charged_higgs_generation/TypeII_joboptions/rivet/routines/pp_dppjj_llvvjj.cc |& tee compile.log
rivet --analysis=pp_dpp_llvvjj --pwd typeII_single_prod_WWjj_100_1.0/Events/run_01/tag_1_pythia8_events.hepmc |& tee run.log
