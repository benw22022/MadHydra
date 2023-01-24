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

def run_gunzip(fpath: str, keep_file: bool=True) -> None:
    """
    Unzips a zipped .gz file (like hepmc.gz produced by pythia)
    args:
        fpath: str - path to file to be unzipped
        keep_file: bool
    """
    
    if os.path.basename(fpath).endswith('.gz'):
        log.info(f"Unzipping: {fpath}")
        unzip_cmd = ['gunzip', fpath]
        if keep_file:
            unzip_cmd.insert(1, '-k')
        launch_process(unzip_cmd)

def run_rivet_routine(config: DictConfig, file_type='*.hepmc*') -> None:    
    
    # Work out input file location
    input_files = glob.glob(os.path.join(config.hydra.runtime.output_dir, config.process.output_dir, 'Events', '*', f'{file_type}'))
    
    print(os.path.join(config.hydra.runtime.output_dir, config.process.output_dir, 'Events', '*', f'{file_type}'))
    print(input_files)
    
    # Decompress files if neccessary
    for fpath in files:
        run_gunzip(fpath)
    
    # Compile routine
    
    # Run routine
    
    # Save result
    