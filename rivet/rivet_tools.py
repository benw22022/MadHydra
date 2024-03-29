"""
Rivet functions
"""

import logger
log = logger.get_logger(__name__)

import os
import glob
from omegaconf import DictConfig
from pathlib import Path
import gzip
import shutil
import inspect
import stat
from source import launch_process, get_files_with_extn
import subprocess

# bash commands for setting up enviroment with rivet
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
        with open(str(fpath).replace(".gz", ""), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return str(fpath).replace('.gz', '')


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

    log.info(f"Compiling {rivet_routine_fpath}")

    # Now compile and run routine, to make life easier write commands to a shell script
    with open("run_rivet.sh", "w") as script:
        script.write("#!/bin/sh\n")
        for cmd in RIVET_SETUP:
            script.write(cmd + "\n")
        script.write(f"rivet-build -r {rivet_routine_fpath} |& tee compile_rivet.log\n")
        script.write(f"rivet --pwd --analysis={routine_name} {hepmc_file} |& tee run_rivet.log\n")

    st = os.stat('run_rivet.sh')
    os.chmod('run_rivet.sh', st.st_mode | stat.S_IEXEC)    
    log.info(f"Running rivet routine {routine_name}")
    subprocess.call(["./run_rivet.sh"], cwd=os.getcwd())


def rivet_analyze_job(config: DictConfig, file_type='.hepmc.gz', routine=None) -> None:    
    """
    Runs rivet analysis for a simulation job
    Args:
        config (DictConfig): hydra config object
        file_type (str, optional): output file extension from MadGraph. Defaults to '*.hepmc.gz'.
        routine (_type_, optional): Name of rivet routine to run, if None then will use routine defined in config. Defaults to None.
    Returns:
        None
    """
    # Work out input file location
    job_dir = config.get("transfer_dir", False)
    if not job_dir or not config.batch:

        # If we are loading a .hydra conf from a previous run
        if config.get("hydra", False):          
            log.info("Reading directory from cfg")
            job_dir = os.path.join(config.hydra.runtime.output_dir, config.process.output_dir)
        
        # Otherwise this is a current MadHydra generation and we're already in the run dir
        else:
            job_dir = os.path.join(os.getcwd(), config.process.output_dir)

    input_files = get_files_with_extn(f"{job_dir}/Events", file_type) # [path for path in Path(job_dir).rglob(f'{file_type}')]

    # If no hepmc files were found
    if len(input_files) < 1:
        log.error(f"No {file_type} files found on path {job_dir}")
        return

    # Override rivet routine in config if requested
    if routine is None:
        routine = config.process.rivet

    # Decompress files if neccessary, compile and run routine
    for fpath in input_files:
        input_fpath = fpath
        if str(fpath).endswith('.gz'):
            try:
                log.info(f"Unzipping file {fpath}")
                input_fpath = run_gunzip(fpath)
            except OSError:
                log.error(f"Could not gunzip file: {fpath}")
                continue
        
        compile_and_run_routine(routine, input_fpath)
    