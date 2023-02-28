"""
Main steering script
"""
import logging
log = logging.getLogger(__name__)

import os
import shutil
import hydra
from hydra.utils import to_absolute_path
from omegaconf import DictConfig, OmegaConf
import source
from rivet import rivet_analyze_job

            
@hydra.main(config_path="config", config_name="config", version_base='1.1')
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
        
        source.launch_process([madgraph_exec, 'proc_card.dat'], 'MadGraph')
        
        if not OmegaConf.is_missing(config, config.process['rivet']):
            log.info(f"Running rivet routine {config.process.rivet}")
            rivet_analyze_job(config)
        
        if config.cleanup:
            log.info("Running cleanup")
            cleanup_cmd = source.get_cleanup_cmd(config) 
            [source.launch_process(cmd.split()) for cmd in cleanup_cmd.split("\n")]
    
if __name__ == "__main__":
    run_generation()