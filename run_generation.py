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

OmegaConf.register_new_resolver("eval", eval) # This will allow us to do arithmetic in our yaml cfgs
            
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
        source.submit_job(config)
    else:
        source.run_local_generation(config)

    log.info("Done")


if __name__ == "__main__":
    run_generation()