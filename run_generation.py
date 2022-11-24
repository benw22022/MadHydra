"""
Run generation of TypeII seesaw model 
https://arxiv.org/abs/1912.08975
"""
import logging
log = logging.getLogger(__name__)

import hydra
from omegaconf import DictConfig
import source
from hydra.utils import get_original_cwd, to_absolute_path
import subprocess
from subprocess import Popen, PIPE, STDOUT
import os
import shutil


def log_subprocess_output(pipe, proc_name=''):
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        log.info(f"{proc_name}: %s", str(line).replace(r"\n", "").strip().replace("b\'", "").rstrip("\'").lstrip("\'"))

@hydra.main(config_path="config", config_name="config")
def run_generation(config: DictConfig) -> None:
    
    if config.batch:
        log.info("Running on batch system")
        source.compose_htc_job(config)
        
        shutil.copyfile(to_absolute_path('config/htc_generation.submit'), os.path.join(os.getcwd(), 'htc_generation.submit'))
        
        print(f"running: {config.process.output_dir}")
        
        os.system(f"condor_submit -batch-name {config.process.output_dir} htc_generation.submit")
        return 

    else:
        log.info("Running local generation")
        log.debug(f"Process is: \n {config.process}")
        source.write_proc_card(config)
        madgraph_exec = os.path.join(to_absolute_path(config.madgraph_dir), 'bin', 'mg5_aMC')
        process = Popen([madgraph_exec, 'proc_card.dat'], stdout=PIPE, stderr=STDOUT)
        with process.stdout:
            log_subprocess_output(process.stdout, 'MadGraph')
        exitcode = process.wait() # 0 means success
        return 
    
if __name__ == "__main__":
    run_generation()