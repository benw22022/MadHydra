#!/usr/bin/env python2
from __future__ import division, print_function

import argparse
import glob
import os

# avoiding athena.py, adapted from https://gitlab.cern.ch/jburr/athena_without_athena 
import AthenaCommon.AtlasUnixStandardJob
from AthenaCommon.AlgSequence import AlgSequence
# from AthenaCommon.AthenaCommonFlags import jobproperties as jps
from AthenaCommon.AppMgr import theApp, ServiceMgr
from AthenaCommon.Debugging import DbgStage
import AthenaCommon.Constants
import AthenaPoolCnvSvc.ReadAthenaPool

from Rivet_i.Rivet_iConf import Rivet_i

from athena_without_athena.setup_argparse import setup_parser

parser = argparse.ArgumentParser("Rivet analysis for the HeavyN search")
parser.add_argument("input", nargs="+", help="input EVNT file(s)")
parser.add_argument("-o", "--output", help='output yoda file')
parser.add_argument("-n", '--nevents', type=int, default=-1, help='process only n events')
parser.add_argument("-l", '--lumi', type=float, default=139, help='luminosity in fb')

setup_parser(parser, default_out=None, msg_level=('--level',), files_input=None)

args = parser.parse_args()

# jps.AthenaCommonFlags.AccessMode = "ClassAccess"
ServiceMgr.MessageSvc.OutputLevel = AthenaCommon.Constants.INFO

theApp.EvtMax = args.nevents

ServiceMgr.EventSelector.InputCollections = args.input

# extract cross section from log.generate
xssum = 0
effsum = 0
nfiles = len(args.input)
for path in args.input:
    inputdir = os.path.split(path)[0]
    with open(os.path.join(inputdir, 'log.generate')) as logfile:
        for line in logfile.read().splitlines():
            if 'MetaData: cross-section (nb)' in line:
                xs = line.split()[-1]
                print("xs = " + str(xs))
                xssum += float(xs) * 1e6 # nb to fb
            elif 'MetaData: GenFiltEff' in line:
                eff = line.split()[-1]
                effsum += float(eff)
            

job = AlgSequence()

rivet = Rivet_i()
rivet.AnalysisPath = os.environ['PWD']
#rivet.Analyses += ['HeavyNRivetAnalysis' ]
rivet.Analyses += ['MLRSM_mumujj']
rivet.RunName = ''
rivet.HistoFile = args.output
rivet.CrossSection = xssum*effsum/nfiles**2 #* args.lumi
print("xssum = " + str(xssum))
print("effsum = " + str(effsum))
print("nfiles = " + str(nfiles))
print("Xsec = " + str(xssum*effsum/nfiles**2 * args.lumi))
rivet.SkipWeights = True
#rivet.IgnoreBeamCheck = True
job += rivet

ServiceMgr.EventSelector.SkipEvents = args.skip_events
# DbgStage.value = 'init'
theApp.run()
print("processed ", ServiceMgr.EventSelector.InputCollections)
theApp.exit()
