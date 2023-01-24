
import os
import argparse

from AthenaCommon.Include import include
from MadGraphControl.MadGraphUtils import *

# General settings
gridpack_mode=False
process_dir=''

stringy = 'aMCPy8EG'

# Merging settings
maxjetflavor=5
ickkw=0
nJetMax=4
ktdurham=30
dparameter=0.4

process = None    # set later explicitly

parameters_runcard = {
    'pdlabel'      : "'lhapdf'",
    'lhaid'        : 260000,
    
}

def run_evgen(runArgs, evgenConfig, opts):
    # scale up number of events so Pythia won't run out of events later on
    
    nevents= runArgs.maxEvents if runArgs.maxEvents > 0 else 5000
    nevents= 1.5 * nevents
    parameters_runcard['nevents']=int(nevents)

    if not is_gen_from_gridpack():
        process_dir = new_process(process)
    else:
        process_dir = MADGRAPH_GRIDPACK_LOCATION
        
    evgenConfig.description = 'test run'
    evgenConfig.keywords+=['BSM','exotic','chargedHiggs'] # list of allowed keywords: https://gitlab.cern.ch/atlas-physics/pmg/infrastructure/mc15joboptions/blob/master/common/evgenkeywords.txt
    evgenConfig.generators += ["aMcAtNlo", "Pythia8", "EvtGen"]
    evgenConfig.process = 'p p ->  same-sign l l j j'
    evgenConfig.inputconfcheck=""
    runArgs.inputGeneratorFile=''.join(['tmp_', stringy, '._00001.events.tar.gz'])

    beamEnergy=-999
    if hasattr(runArgs,'ecmEnergy'):
        beamEnergy = runArgs.ecmEnergy / 2.
    else: 
        raise RuntimeError("No center of mass energy found.")
        
    modify_run_card(
        run_card_input=get_default_runcard(process_dir),
        process_dir=process_dir,
        runArgs=runArgs,
        settings=parameters_runcard
    )
    
    print(parameters_paramcard)

    modify_param_card(
        param_card_input='/'.join([process_dir, 'Cards', 'param_card.dat']),
        process_dir=process_dir,
        params={}#parameters_paramcard
    )
                    
    print_cards()

    generate(
        process_dir=process_dir,
        grid_pack=gridpack_mode,
        runArgs=runArgs
    )
    
    arrange_output(
        process_dir=process_dir,
        runArgs=runArgs,
        lhe_version=3,
        saveProcDir=True
    )
    
    # back to single core
    check_reset_proc_number(opts)
    
    include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
    include("Pythia8_i/Pythia8_aMcAtNlo.py")  
    
    # Execute event gen
    run_evgen(runArgs, evgenConfig, opts)

# gen_command =f"Gen_tf.py --ignorePatterns=\"attempt to add a duplicate\" --ecmEnergy=13000. --maxEvents={str(args.nevents)} --firstEvent=1 --randomSeed={str(args.randomN)} --outputEVNTFile={outfile} --jobConfig=${PWD}/JO --asetup=''"