from MadGraphControl.MadGraphUtils import *
     
import os

evgenConfig.generators += ["MadGraph", "EvtGen"]
evgenConfig.keywords = ['SM', 'diboson', 'VBS', 'WW', 'sameSign', 'electroweak', '2lepton', '2jet']
evgenConfig.contact = ['Karolos Potamianos <karolos.potamianos@cern.ch>', 'Cyril Becot <cyril.becot@cern.ch>', 'Oleg Kuprash <oleg.kuprash@cern.ch>', 'Carsten Bittrich <carsten.bittrich@cern.ch>']

evgenConfig.nEventsPerJob = 10000

safetyFactor = 1.1
nevents = int(runArgs.maxEvents*safetyFactor if runArgs.maxEvents>0 else safetyFactor*evgenConfig.nEventsPerJob)

required_accuracy = 0.01

# Todo for new paradigm: switch to groups in naming, i.e. not based on DSIDs

"""
Handling JOs for 12 DSIDs: 500986-500994
500986: EW6  LO Py8
500987: EW4  LO Py8
500988: INT  LO Py8
500989: EW6  LO  H7
500990: EW4  LO  H7
500991: INT  LO  H7
500992: EW6 NLO  H7
500993: EW4 NLO  H7
500994: INT NLO  H7
"""

isEW6 = [ 500986, 500989, 500992 ]
isEW4 = [ 500987, 500990, 500993 ]
isINT = [ 500988, 500991, 500994 ]

isPythia = [ 500986, 500987, 500988 ]
isHerwig = [ 500989, 500990, 500991, 500992, 500993, 500994 ]

isLO  = [ 500986, 500987, 500988, 500989, 500990, 500991 ]
isNLO = [ 500992, 500993, 500994 ]

# Use dipole shower for LO and Pythia
useDipoleRecoil = list(set(isLO).intersection(isPythia))
# Use dipole shower for LO and Herwig
useDipoleShower = list(set(isLO).intersection(isHerwig))

#thisDSID = int(runArgs.jobConfig[0])
thisDSID = int(os.path.basename(runArgs.jobConfig[0]))

if thisDSID in isNLO and thisDSID in (useDipoleRecoil+useDipoleShower):
  errMsg = "Cannot use dipole recoil or dipole shower with NLO"
  print errMsg
  raise RuntimeError(errMsg)

if thisDSID in isPythia: evgenConfig.generators += ["Pythia8"]
elif thisDSID in isHerwig: evgenConfig.generators += ["Herwig7"]

gridpack_mode=False

# Processes


if thisDSID in isEW6:
  if thisDSID in isLO:
    runName = 'lvlvjj_ss_EW6'
    description = 'MadGraph_' + runName
    process = """
  import model sm-no_b_mass
  define l+ = e+ mu+ ta+
  define vl = ve vm vt
  define l- = e- mu- ta-
  define vl~ = ve~ vm~ vt~
  define p = g u c d s u~ c~ d~ s~ b b~
  define j = g u c d s u~ c~ d~ s~ b b~
  generate p p > l+ vl l+ vl j j QED==6 QCD=0 @1
  add process p p > l- vl~ l- vl~ j j QED==6 QCD=0 @1
  output -f"""
  elif thisDSID in isNLO:
    runName = 'lvlvjj_ss_EW6_NLO'
    description = 'MadGraph_' + runName
    process = """
  set complex_mass_scheme
  import model loop_qcd_qed_sm_Gmu-atlas
  define p = g u c d s u~ c~ d~ s~ b b~
  define j = g u c d s u~ c~ d~ s~ b b~
  generate p p > w+ w+  j j QED=4 QCD=0 [QCD] @1
  add process p p > w- w-  j j QED=4 QCD=0 [QCD] @1
  output -f"""

elif thisDSID in isINT:
  if thisDSID in isLO:
    runName = 'lvlvjj_ss_INT'
    description = 'MadGraph_' + runName
    process = """
  import model sm-no_b_mass
  define l+ = e+ mu+ ta+
  define vl = ve vm vt
  define l- = e- mu- ta-
  define vl~ = ve~ vm~ vt~
  define p = g u c d s u~ c~ d~ s~ b b~
  define j = g u c d s u~ c~ d~ s~ b b~
  generate p p > l+ vl l+ vl j j QCD^2==2 @1
  add process p p > l- vl~ l- vl~ j j QCD^2==2 @1
  output -f"""
  elif thisDSID in isNLO:
    runName = 'lvlvjj_ss_INT_NLO'
    description = 'MadGraph_' + runName
    process = """
  set complex_mass_scheme
  import model loop_qcd_qed_sm_Gmu-atlas
  define p = g u c d s u~ c~ d~ s~ b b~
  define j = g u c d s u~ c~ d~ s~ b b~
  generate p p > w+ w+ j j QCD^2==2 [QCD] @1
  add process p p > w- w- j j QCD^2==2 [QCD] @1
  output -f"""

elif thisDSID in isEW4:
  if thisDSID in isLO:
    runName = 'lvlvjj_ss_EW4'
    description = 'MadGraph_' + runName
    process = """
  import model sm-no_b_mass
  define l+ = e+ mu+ ta+
  define vl = ve vm vt
  define l- = e- mu- ta-
  define vl~ = ve~ vm~ vt~
  define p = g u c d s u~ c~ d~ s~ b b~
  define j = g u c d s u~ c~ d~ s~ b b~
  generate p p > l+ vl l+ vl j j QED=4 QCD=2 @1
  add process p p > l- vl~ l- vl~ j j QED=4 QCD=2 @1
  output -f"""
  elif thisDSID in isNLO:
    runName = 'lvlvjj_ss_EW4_NLO'
    description = 'MadGraph_' + runName
    process = """
  set complex_mass_scheme
  import model loop_qcd_qed_sm_Gmu-atlas
  define p = g u c d s u~ c~ d~ s~ b b~
  define j = g u c d s u~ c~ d~ s~ b b~
  generate p p > w+ w+ j j QED=2 QCD=2 [QCD] @1
  add process p p > w- w- j j QED=2 QCD=2 [QCD] @1
  output -f"""

else:
  raise RuntimeError("runNumber %i not recognised in these jobOptions." % thisDSID)

evgenConfig.description = description

#Fetch default NLO run_card.dat and set parameters
extras = { 'lhe_version'  :'3.0',
           'pdlabel'      :"'lhapdf'",
           'lhaid'        : 260000,
           'bwcutoff'     :'1000',
           'nevents'      :nevents,
           'dynamical_scale_choice': 0,
           'maxjetflavor': 5, 
           'ptl': 4.0, 
           'ptj': 15.0, 
           'etal': 3.0, 
           'etaj': 5.5, 
           'drll': 0.2, 
           'systematics_program': 'systematics',  
           'systematics_arguments': "['--mur=0.5,1,2', '--muf=0.5,1,2', '--dyn=-1,1,2,3,4', '--weight_info=MUR%(mur).1f_MUF%(muf).1f_PDF%(pdf)i_DYNSCALE%(dyn)i', '--pdf=errorset,NNPDF30_nlo_as_0119@0,NNPDF30_nlo_as_0117@0,CT14nlo@0,MMHT2014nlo68clas118@0,PDF4LHC15_nlo_30_pdfas']"}

process_dir = new_process(process)

if not is_gen_from_gridpack(): # When generating the gridpack
  # Todo: replace in line to be safer
  os.system("cp collect_events.f "+process_dir+"/SubProcesses")


if thisDSID in isLO:
  if not is_gen_from_gridpack(): # When generating the gridpack
    os.system("cp setscales_lo.f  "+process_dir+"/SubProcesses/setscales.f")
  lo_extras = { 'asrwgtflavor': 5, 
                'auto_ptj_mjj': False, 
                'cut_decays': True, 
                'ptb': 15.0, 
                'etab': 5.5, 
                'dral': 0.1, 
                'drbb': 0.2, 
                'drbj': 0.2, 
                'drbl': 0.2, 
                'drjj': 0.2,
                'drjl': 0.2, 
                'mmll': 0.0, 
                'gridpack': '.true.',
                'use_syst': True, 
              }
  extras.update(lo_extras)

elif thisDSID in isNLO:
  if not is_gen_from_gridpack(): # When generating the gridpack
    # See https://answers.launchpad.net/mg5amcnlo/+faq/2720
    # Todo: replace in line to be safer
    os.system("cp FKS_params.dat "+process_dir+"/Cards/")
    os.system("cp setscales_nlo.f  "+process_dir+"/SubProcesses/setscales.f")

  nlo_extras = { 'lhe_version'   :'3.0',
           'pdlabel'      : "'lhapdf'",
           'lhaid'         : 260000,
           'parton_shower' : 'PYTHIA8' if thisDSID in isPythia else 'HERWIGPP',
           'jetalgo'       : -1,
           'jetradius'     : 0.4,
           #'ptj'           : 20.0,
           #'etaj'          : -1,
           'ickkw'         : 0,
           'store_rwgt_info' : '.true.',
           'event_norm' : 'sum',
           'reweight_scale': '.true.',
           'reweight_PDF'  : '.true.',
           'dynamical_scale_choice': 10  
         }
  extras.update(nlo_extras)

  required_accuracy = 0.001

  if not is_gen_from_gridpack(): # When generating the gridpack
    madspin_card = process_dir+"/Cards/madspin_card.dat"
    # If one wants to use runName one needs to change MADGRAPH_RUN_NAME anyhow (tbc)
    ms_import = process_dir+"/Events/"+MADGRAPH_RUN_NAME+"/events.lhe"
    ms_dir = "MadSpin"
  else:
    gridpack_dir = "madevent/"
    madspin_card = gridpack_dir+"Cards/madspin_card.dat"
    # If one wants to use runName one needs to change MADGRAPH_RUN_NAME anyhow (tbc)
    ms_import = gridpack_dir+"Events/"+MADGRAPH_RUN_NAME+"/events.lhe"
    ms_dir = gridpack_dir+"MadSpin"

  if os.access(madspin_card,os.R_OK):
    os.unlink(madspin_card)
  mscard = open(madspin_card,'w')  
  mscard.write("""#************************************************************
#*                        MadSpin                           *
#*                                                          *
#*    P. Artoisenet, R. Frederix, R. Rietkerk, O. Mattelaer *
#*                                                          *
#*    Part of the MadGraph5_aMC@NLO Framework:              *
#*    The MadGraph5_aMC@NLO Development Team - Find us at   *
#*    https://server06.fynu.ucl.ac.be/projects/madgraph     *
#*                                                          *
#************************************************************
#Some options (uncomment to apply)
#
# set seed 1
# set Nevents_for_max_weigth 75 # number of events for the estimate of the max. weight
# set BW_cut 15                # cut on how far the particle can be off-shell
 set max_weight_ps_point 400  # number of PS to estimate the maximum for each event
#
import %s
set ms_dir %s
set seed %i
# specify the decay for the final state particles
define l+ = e+ mu+ ta+
define l- = e- mu- ta-
define v = ve vm vt
define v~ = ve~ vm~ vt~
decay w+ > l+ v                                                                   
decay w- > l- v~          
# running the actual code
launch"""%(ms_import, ms_dir, 10000000+int(runArgs.randomSeed)))
  mscard.close()

modify_run_card(runArgs=runArgs,process_dir=process_dir,settings=extras)
#modify_run_card(process_dir=process_dir,runArgs=runArgs,settings=settings)

generate(process_dir=process_dir,grid_pack=gridpack_mode,runArgs=runArgs,required_accuracy=required_accuracy)
arrange_output(runArgs=runArgs,process_dir=process_dir,lhe_version=3,saveProcDir=True)

#### Shower                    

if thisDSID in isPythia:
  include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
  include("Pythia8_i/Pythia8_aMcAtNlo.py")
  # fix for VBS processes (requires version>=8.230)
  genSeq.Pythia8.Commands += [ "SpaceShower:pTmaxMatch=1",
                               "SpaceShower:pTmaxFudge=1",
                               "SpaceShower:MEcorrections=off",
                               "TimeShower:pTmaxMatch=1",
                               "TimeShower:pTmaxFudge=1",
                               "TimeShower:MEcorrections=off",
                               "TimeShower:globalRecoil=on",
                               "TimeShower:limitPTmaxGlobal=on",
                               "TimeShower:nMaxGlobalRecoil=1",
                               "TimeShower:globalRecoilMode=2",
                               "TimeShower:nMaxGlobalBranch=1",
                               "TimeShower:weightGluonToQuark=1",
                               "Check:epTolErr=1e-2"]

  if thisDSID in useDipoleRecoil:
    genSeq.Pythia8.Commands += [ 'SpaceShower:dipoleRecoil = on' ]


  include("Pythia8_i/Pythia8_MadGraph.py")
  # include("Pythia8_i/Pythia8_ShowerWeights.py")

elif thisDSID in isHerwig:
  from Herwig7_i.Herwig7_iConf import Herwig7
  from Herwig7_i.Herwig72ConfigLHEF import Hw72ConfigLHEF

  genSeq += Herwig7()
  Herwig7Config = Hw72ConfigLHEF(genSeq, runArgs)

# configure Herwig7
  Herwig7Config.me_pdf_commands(order="NLO", name="NNPDF30_nlo_as_0118")
  Herwig7Config.tune_commands()
  lhe_filename = runArgs.inputGeneratorFile.split('.tar.gz')[0]+'.events'
  Herwig7Config.lhef_mg5amc_commands(lhe_filename=lhe_filename, me_pdf_order="NLO")
  
  if thisDSID in useDipoleShower:
     dipoleShowerCommands = """
     cd /Herwig/EventHandlers
     set EventHandler:CascadeHandler /Herwig/DipoleShower/DipoleShowerHandler
     cd /Herwig/DipoleShower
     do DipoleShowerHandler:AddVariation isr:muRfac=2.0_fsr:muRfac=2.0 2.0 2.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=2.0_fsr:muRfac=1.0 2.0 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=2.0_fsr:muRfac=0.5 2.0 0.5 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=2.0" 1.0 2.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=0.5" 1.0 0.5 All
     do DipoleShowerHandler:AddVariation isr:muRfac=0.5_fsr:muRfac=2.0" 0.5 2.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=0.5_fsr:muRfac=1.0 0.5 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=0.5_fsr:muRfac=0.5 0.5 0.5 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.75_fsr:muRfac=1.0 1.75 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.5_fsr:muRfac=1.0 1.5 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.25_fsr:muRfac=1.0 1.25 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=0.625_fsr:muRfac=1.0" 0.625 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=0.75_fsr:muRfac=1.0 0.75 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=0.875_fsr:muRfac=1.0 0.875 1.0 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=1.75 1.0 1.75 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=1.5 1.0 1.5 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=1.25 1.0 1.25 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=0.625 1.0 0.625 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=0.75 1.0 0.75 All
     do DipoleShowerHandler:AddVariation isr:muRfac=1.0_fsr:muRfac=0.875 1.0 0.85 All
     """
     print dipoleShowerCommands
     Herwig7Config.add_commands(dipoleShowerCommands)
  
# add EvtGen
  include("Herwig7_i/Herwig71_EvtGen.py")

# run Herwig7
  Herwig7Config.run()
