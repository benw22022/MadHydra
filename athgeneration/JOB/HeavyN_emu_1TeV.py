import subprocess
retcode = subprocess.call(['get_files', '-jo', 'HeavyNCommon.py'])
if retcode != 0:
    raise IOError('could not locate HeavyNCommon.py')

# import HeavyNCommon as HeavyNCommon

import os

from AthenaCommon.Include import include
from MadGraphControl.MadGraphUtils import *

# General settings
gridpack_mode=True
process_dir=''

stringy = 'aMCPy8EG_NNPDF3NLO_HeavyN_VBS_NLO'

# Merging settings
maxjetflavor=5
ickkw=0
nJetMax=4
ktdurham=30
dparameter=0.4

process = None                  # set later explicitly

available_processes = {
    'pythia':'p p > mu+- mu+- j j',   # for pythia
    'mumuchannel':'''
import model SM_HeavyN_NLO

define p = u c d s u~ c~ d~ s~ g
define j = p
define mu = mu+ mu-

generate p p > mu+ mu+ j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]
add process p p > mu- mu- j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]

output -f
    ''',
    'eechannel':'''
import model SM_HeavyN_NLO

define p = u c d s u~ c~ d~ s~ g
define j = p
define e = e+ e-

generate p p > e+ e+ j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]
add process p p > e- e- j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]

output -f
    ''',
    'tautauchannel':'''
import model SM_HeavyN_NLO

define p = u c d s u~ c~ d~ s~ g
define j = p
define ta = ta+ ta-

generate p p > ta+ ta+ j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]
add process p p > ta- ta- j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]

output -f
    ''',

'emuchannel':'''
import model SM_HeavyN_NLO

define p = u c d s u~ c~ d~ s~ g
define j = p

generate p p > e+ mu+ j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]
add process p p > e- mu- j j $$ n1 n2 n3 w+ w- QED=4 QCD=0 [QCD]

output -f
    '''

}


parameters_runcard = {
    'pdlabel'      : "'lhapdf'",
    'lhaid'        : 260000,
    'jetalgo':-1,
    'jetradius':0.4,
    'etaj':5.2,                 # generate until slightly out of detector acceptance
    'etal':2.9,                 # maybe some smearing effects could bring them into acceptance
    'mll':15,
    'ptj':15,
    'ptl':10,
    'parton_shower':'PYTHIA8'
}

parameters_paramcard = {
    'mass':{
        'mN1':1e10,
        'mN2':1e10,
        'mN3':1e10,
    },
    'numixing':{
        'VeN1':0,
        'VeN2':0,
        'VeN3':0,
        'VmuN1':0,
        'VmuN2':0,
        'VmuN3':0,
        'VtaN1':0,
        'VtaN2':0,
        'VtaN3':0,
    }
}

# modify_run_card(
#     process_dir=process_dir,
#     runArgs=runArgs,
#     settings=extras_runcard
# )
# modify_config_card(
#     process_dir=process_dir
# )

# generate(
#     process_dir=process_dir,
#     grid_pack=gridpack_mode,
#     runArgs=runArgs
#     # nevents=nevents
# )


def run_evgen(runArgs, evgenConfig, opts):
    # scale up number of events so Pythia won't run out of events later on
    nevents=runArgs.maxEvents if runArgs.maxEvents > 0 else 5000
    nevents=1.5*nevents
    parameters_runcard['nevents']=int(nevents)

    process = available_processes['emuchannel']
    parameters_paramcard['mass']['mN1'] = 1000
    parameters_paramcard['numixing']['VmuN1'] = 1
    parameters_paramcard['numixing']['VeN1'] = 1


    if not is_gen_from_gridpack():
        process_dir = new_process(process)
    else:
        process_dir = MADGRAPH_GRIDPACK_LOCATION
        
        
    evgenConfig.description = 'test run'
    evgenConfig.keywords+=['BSM','exotic','neutrino', 'VBS'] # list of allowed keywords: https://gitlab.cern.ch/atlas-physics/pmg/infrastructure/mc15joboptions/blob/master/common/evgenkeywords.txt
    evgenConfig.generators += ["aMcAtNlo", "Pythia8", "EvtGen"]
    evgenConfig.process = 'pp -> mu mu j j'
    evgenConfig.inputconfcheck=""
    evgenConfig.contact = ["Jonas Neundorf <jonas.neundorf@desy.de>"]
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
    
    modify_param_card(
        param_card_input='/'.join([process_dir, 'Cards', 'param_card.dat']),
        process_dir=process_dir,
        params=parameters_paramcard
    )
    
    # modify FKS card to ignore IRPoleCheckThreshold
    print 'Modifying SKS Param Card, printing old stuff'
    with open('/'.join([process_dir, 'Cards', 'FKS_params.dat']), 'r+') as fks:
        retlines = []
        for line in fks:
            if '1.0d-5' in line and not line.startswith('!'): # at the time of writing, only param with this value
                retlines.append('-1.0d0\n')
            else:
                retlines.append(line)

        fks.seek(0)                 # set "cursor" at beginning of file
        fks.write(''.join(retlines))
        fks.truncate()              # and cut off the end
                
                
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

    run_evgen(runArgs, evgenConfig, opts)


# HeavyNCommon.process = HeavyNCommon.available_processes['emuchannel']
# HeavyNCommon.parameters_paramcard['mass']['mN1'] = 1000
# HeavyNCommon.parameters_paramcard['numixing']['VmuN1'] = 1
# HeavyNCommon.parameters_paramcard['numixing']['VeN1'] = 1



# HeavyNCommon.run_evgen(runArgs, evgenConfig, opts)
