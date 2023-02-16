"""
Rivet functions
"""

import logging
log = logging.getLogger(__name__)

import os
import glob
from omegaconf import DictConfig
from source import launch_process
from typing import List
import gzip
import shutil
import inspect
import stat


RIVET_SETUP = (r'export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
               r'source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh',
               r'asetup 23.6.1, AthGeneration',
               r'lsetup "views LCG_102b_ATLAS_2 x86_64-centos7-gcc11-opt"')

def run_gunzip(fpath: str) -> None:
    """
    Unzips a zipped .gz file (like hepmc.gz produced by pythia)
    args:
        fpath: str - path to file to be unzipped
        keep_file: bool - keep compressed file after unzipping
    returns:
        None
    """
    
    with gzip.open(fpath, 'rb') as f_in:
        with open(fpath.replace('.gz', ''), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return fpath.replace('.gz', '')

def compile_and_run_routine(routine_name: str, hepmc_file: str) -> None:
    """
    Script to compile and run rivet routine
    Requires cvmfs and Athena 23.6.1
    args:
        routine_name: str - name of rivet routine to run; <routine>.cpp must be in rivet/routines
        hepmc_file: str - path to the hepmc to be analyzed with this rivet routine
    returns:
        None
    """

    # Locate module directory (need to do this to reliably locate routines)
    this_filename = inspect.getframeinfo(inspect.currentframe()).filename
    this_filepath = os.path.dirname(os.path.abspath(this_filename))

    # Grab routine name from config and extract fpath to C++ file
    rivet_routine_fpath = os.path.join(this_filepath, "routines", f"{routine_name}.cc")

    # Now compile and run routine, to make life easier write commands to a shell script
    with open("run_rivet.sh", "w") as script:
        script.write("#!/bin/sh\n")
        for cmd in RIVET_SETUP:
            script.write(cmd + "\n")
        script.write(f"rivet-build -r {rivet_routine_fpath} |& tee compile.log\n")
        script.write(f"rivet --analysis={routine_name} --pwd {hepmc_file} |& tee run.log\n")

    st = os.stat('run_rivet.sh')
    os.chmod('run_rivet.sh', st.st_mode | stat.S_IEXEC)    
    # launch_process(["./run_rivet.sh"], "Rivet")
    os.system()


def rivet_analyze_job(config: DictConfig, file_type='*.hepmc.gz', keep_gz_files=False, routine=None) -> None:    

    # Work out input file location
    job_dir = config.get("transfer_dir", False)
    if not job_dir:
        job_dir = config.hydra.runtime.output_dir
    else:
        job_dir = os.path.join(job_dir, config.process.output_dir)

    events_dir = os.path.join(job_dir, config.process.output_dir, 'Events')
    input_files = glob.glob(os.path.join(events_dir, '*', f'{file_type}'))

    # If no hepmc files were found
    if len(input_files) < 1:
        log.error(f"No hepmc files found on path {events_dir}")
        return

    # Override rivet routine in config if requested
    if routine is None:
        routine = config.process.rivet

    # Decompress files if neccessary, compile and run routine
    for fpath in input_files:
        input_fpath = fpath
        if fpath.endswith('.gz'):
            try:
                log.info(f"Unzipping file {fpath}")
                input_fpath = run_gunzip(fpath)
            except OSError:
                log.error(f"Could not gunzip file: {fpath}")
                continue

        compile_and_run_routine(routine, input_fpath)
    