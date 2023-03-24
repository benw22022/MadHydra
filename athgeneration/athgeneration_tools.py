"""
Tools for running using AthGeneration
"""
import logger
log = logger.get_logger(__name__)

import os
from omegaconf import DictConfig
import inspect
import subprocess

# bash commands for setting up enviroment with AthGeneration
ATHGEN_SETUP = (r'export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase',
                r'source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh',
                r'asetup AthGeneration, 22.6.24',
                r'export ATHENA_PROC_NUMBER=12')


def generate_dsid(config: DictConfig) -> None:

    # Locate module directory (need to do this to reliably locate routines)
    this_filename = inspect.getframeinfo(inspect.currentframe()).filename
    this_filepath = os.path.dirname(os.path.abspath(this_filename))

    # Grab routine name from config and extract fpath to C++ file
    script_path = os.path.join(this_filepath, "gen.sh")

    # Make sure that the user has actually specified a dsid in their config
    try:
        dsid = config.process.dsid
    except KeyError:
        log.error("Error: could not find dsid in process config file")

    log.info(f"Generating sample of {config.parameters.nevents} for {config.process.process_name}:{config.process.dsid}")
    
    # If we're running in debug mode we can stop here
    if config.debug: return

    # Launch script which calls Gen_tf.py
    subprocess.call([f"{script_path}", f"{dsid}", f"{config.process.output_dir}", f"{config.parameters.nevents}"], cwd=os.getcwd())
