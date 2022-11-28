import logging
log = logging.getLogger(__name__)

import os
from omegaconf import DictConfig
from hydra.utils import get_original_cwd, to_absolute_path

def write_proc_card(config : DictConfig) -> None:
    """
    Writes the proc card for the process given in the config file
    args: 
        config: DictConfig - config object
    returns:
        None
    """
    
    with open('proc_card.dat', 'w') as proc_card:
        for line in config.process.gen_cmd:
            proc_card.write(f'{line}\n')
            
def compose_htc_job(config: DictConfig) -> None:
    """
    Makes scripts to run job on batch system
    Yes this does look like an awful way of writing files, seems that there is a bug in how htc processes
    bash scripts. Writing out the lines one by one, or writing a block comment to file causes htc to not execute
    the contents of the script at all. You have to write the entire script as a single python string.
    args: 
        config: DictConfig - config object
    returns:
        None
    """
    
    
    write_proc_card(config)
    madgraph_exec = os.path.join(to_absolute_path(config.madgraph_dir), 'bin', 'mg5_aMC')
    cwd = os.getcwd()
    
    # Create python env for MadGraph
    with open('setup_and_run.sh', 'w') as run_script:
        cmd = f"#!/bin/bash \n \
source cd {cwd} \n \
chmod +x  run_generation.sh \n \
source ./run_generation.sh" 
        run_script.write(cmd) # yes | conda create -n python39 python=3.9 \n \
    
    # Execute MadGraph 
    with open('run_generation.sh', 'w') as run_script:
        cmd = f"#!/bin/bash \n \
eval \"$(conda shell.bash hook)\" \n \
conda activate python39 \n 
python3 {madgraph_exec} proc_card.dat | tee log.generate"
        run_script.write(cmd) # pip3 install six --user \n \
                
