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
import os



@hydra.main(config_path="config", config_name="config")
def run_generation(config: DictConfig) -> None:
    
    if config.batch:
        # TODO add batch running
        log.info("Running on batch system")
        return 

    else:
        log.info("Running local generation")
        log.debug(f"Process is: \n {config.process}")
        source.write_proc_card(config)
        madgraph_exec = os.path.join(to_absolute_path(config.madgraph_dir), 'bin', 'mg5_aMC')
        subprocess.run([madgraph_exec, 'proc_card.dat'])
    
if __name__ == "__main__":
    run_generation()