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
    
    # Create script to execute MadGraph 
    with open('run_generation.sh', 'w') as run_script:
        cmd = f"#!/bin/bash \n \
eval \"$(conda shell.bash hook)\" \n \
conda activate python39 \n \
python3 {madgraph_exec} {cwd}/proc_card.dat | tee log.generate \n"
        
        # Run clean up of MG dir (avoid running into disk quota limits!)
        if config.clean_up:
            cmd += f"rm -r tmp* \n"
            cmd += f"rm -r py.py \n"
            cmd += f"rm -r {config.process.output_dir}/bin \n"
            cmd += f"rm -r {config.process.output_dir}/Source \n"
            cmd += f"rm -r {config.process.output_dir}/lib \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*.f \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*.inc \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/Makefile \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/done \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/.txt \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/.mg \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/.sh \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/.dat \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/randinit \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/proc_characteristics \n"
            
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/*.f \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/*.inc \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/*.sym \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/Makefile \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/.txt \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/.mg \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/.sh \n"
            cmd += f"rm -r {config.process.output_dir}/SubProcesses/*/.dat \n"
         
        # Transfer files back to launch dir or transfer them to another location like eos   
        if config.transfer_files:
            cmd += f"\nmv {config.process.output_dir} {config.transfer_dir}"
        else:
            cmd += f"\n mv * {cwd}"
        
        # Write script 
        run_script.write(cmd)