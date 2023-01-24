"""
Main steering script
"""
import logging
log = logging.getLogger(__name__)

import os
import shutil
import subprocess
from subprocess import Popen, PIPE, STDOUT
import hydra
from hydra.utils import get_original_cwd, to_absolute_path
from omegaconf import DictConfig
import source
from typing import List
import re

def escape_ansi(line):
    ansi_escape =re.compile(r'(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]')
    return ansi_escape.sub('', line)

def log_subprocess_output(pipe, proc_name: str='') -> None:
    """
    Logs output from subprocess to log file
    args:
        pipe: stdout buffer to read from
        proc_name: str (default='')- a string to optionally label the process with
    returns:
        None
    """
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        cleaned_line = escape_ansi(line)
        cleaned_line = str(cleaned_line).replace(r"\n", "").strip().replace("b\"", "").rstrip("\'").lstrip("\'")
        log.info(f"{proc_name}: %s", cleaned_line)

def launch_process(cmd_list: List[str], proc_name: str='') -> int:
    """
    Launches a subprocess, just a little function to avoid some of the boilerplate
    args:
        cmd_list: List[str] - a list of commands/options to execute. e.g "rm -r aDir" would be ['rm', '-r', 'aDir']
        proc_name: str (default='')- a string to optionally label the process with when logging
    returns:
        int: exit code (0 means success)
    """
    
    if len(cmd_list) == 0:
        return 0
    
    process = Popen(cmd_list, stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess_output(process.stdout, proc_name)
    return process.wait()
            
@hydra.main(config_path="config", config_name="config")
def run_generation(config: DictConfig) -> None:
    """
    Main steering script and hydra entry point
    args:
        config: DictConfig - Hydra configuration, automagically parsed by @hydra.main decorator
    returns:
        None
    """
    
    if config.batch:
        log.info("Running on batch system")
        source.compose_htc_job(config)
        
        shutil.copyfile(to_absolute_path('config/htc_generation.submit'), os.path.join(os.getcwd(), 'htc_generation.submit'))
        
        log.info(f"running: {config.process.output_dir}")
        
        if config.debug:
            return
        
        os.system(f"condor_submit -batch-name {config.process.output_dir} htc_generation.submit")
        return 

    else:
        log.info("Running local generation")
        log.debug(f"Process is: \n {config.process}")
        source.write_proc_card(config)
        madgraph_exec = os.path.join(to_absolute_path(config.madgraph_dir), 'bin', 'mg5_aMC')
        
        if config.debug:
            return
        
        launch_process([madgraph_exec, 'proc_card.dat'], 'MadGraph')
        
        if config.cleanup:
            log.info("Running cleanup")
            cleanup_cmd = source.get_cleanup_cmd(config) 
            os.system("pwd")
            print([cmd.split() for cmd in cleanup_cmd.split("\n")])
            
            [launch_process(cmd.split()) for cmd in cleanup_cmd.split("\n")]
    
if __name__ == "__main__":
    run_generation()