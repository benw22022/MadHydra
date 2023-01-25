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

RIVET_SETUP = ('setupATLAS', 'lsetup root', 'asetup 21.6.1, AthGeneration','source ${LCG_RELEASE_BASE}/LCG_88/MCGenerators/rivet/3.1.2/${LCG_PLATFORM}/rivetenv.sh')

def run_gunzip(fpath: str, keep_file: bool=True) -> None:
    """
    Unzips a zipped .gz file (like hepmc.gz produced by pythia)
    args:
        fpath: str - path to file to be unzipped
        keep_file: bool - keep compressed file after unzipping
    returns:
        None
    """
    
    if os.path.basename(fpath).endswith('.gz'):
        log.info(f"Unzipping: {fpath}")
        unzip_cmd = ['gunzip', fpath]
        if keep_file:
            unzip_cmd.insert(1, '-c')
            unzip_cmd += ['>', fpath.replace('.gz', "")]

        # TODO: Make this work without shell=True
        launch_process(" ".join(unzip_cmd), shell=True)

def compile_routine(config: DictConfig) -> None:

    # TODO: Bug here - will not work for tests at different locations
    routines_dir = os.path.join(config.hydra.runtime.cwd, 'rivet', 'routines')
    routine = os.path.join(routines_dir, f'{config.process.rivet}.cc')
    compile_script = f"compile_rivet_routine_{config.process.rivet}.sh"

    with open(compile_script, 'w') as file:
        file.write("echo hello\n")
        for cmd in RIVET_SETUP:
            file.write(f"{cmd}\n")
        file.write(f"rivet-build {routine}\n")
        
    xperm_cmd = ["chmod", "+x", compile_script]
    compile_cmd = f"./{compile_script}"

    # retcode = launch_process(compile_cmd)
    retcode = launch_process(" ".join(xperm_cmd), shell=True)
    retcode = launch_process(" ".join(compile_cmd), shell=True)
    os.system(compile_cmd)

    # if retcode != 0:
    #     log.error(f"Compilation error! Rivet routine {routine} failed to compile")
    #     raise RuntimeError

def run_rivet_routine(config: DictConfig, file_type='*.hepmc*', keep_gz_files=False) -> None:    
    
    # Work out input file location
    input_files = glob.glob(os.path.join(config.hydra.runtime.output_dir, config.process.output_dir, 'Events', '*', f'{file_type}'))

    print(input_files)    
    
    # Decompress files if neccessary
    for fpath in input_files:
        print("Unzipping ", fpath)
        # run_gunzip(fpath, keep_file=keep_gz_files)

    # Compile routine
    compile_routine(config)
    
    

    # Run routine
    
    # Save result
    