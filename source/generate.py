import logger
log = logger.get_logger(__name__)

import os
import shutil
import inspect
import glob
import hydra
from omegaconf import DictConfig, OmegaConf
import distutils
from distutils.dir_util import copy_tree
from hydra.utils import get_original_cwd, to_absolute_path

import source.utils as utils
import rivet



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


def run_cleanup(config: DictConfig) -> None:
    """
    Deletes unncessary files that cluster up output directory 
    Need to prevent file system quoata filling up
    args:
        config: DictConfig - config object
    returns:
        None
    """

    utils.delete_files_with_extn(config.process.output_dir, '.f')
    utils.delete_files_with_extn(config.process.output_dir, '.o')
    utils.delete_files_with_extn(config.process.output_dir, '.inc')
    utils.delete_files_with_extn(config.process.output_dir, '.txt')
    utils.delete_files_with_extn(config.process.output_dir, '.mg')
    utils.delete_files_with_extn(config.process.output_dir, '.sh')
    utils.delete_files_with_extn(config.process.output_dir, '.dat')
    shutil.rmtree(f"{config.process.output_dir}/bin")
    shutil.rmtree(f"{config.process.output_dir}/Source")
    shutil.rmtree(f"{config.process.output_dir}/lib")


def run_local_generation(config: DictConfig) -> None:
    """
    Runs madgraph generation in current working directory
    args:
        config: DictConfig - config object
    returns:
        None
    """

    log.info("Running local generation")
    log.debug(f"Process is: \n {config.process}")

    # Some hacky stuff to workout where MadGraph is 
    this_filename = inspect.getframeinfo(inspect.currentframe()).filename
    this_filepath = os.path.dirname(os.path.dirname(os.path.abspath(this_filename)))
    mg_path = os.path.join(this_filepath, config.madgraph_dir)
    madgraph_exec = os.path.join(mg_path, 'bin', 'mg5_aMC')

    write_proc_card(config)
    
    if config.debug:
        return
    
    utils.launch_process([madgraph_exec, 'proc_card.dat'], 'MadGraph', logfile='log.generate')
    
    if not OmegaConf.is_missing(config, config.process['rivet']):
        log.info(f"Running rivet routine {config.process.rivet}")
        rivet.rivet_analyze_job(config)
    
    if config.cleanup:
        log.info("Running cleanup")
        run_cleanup(config)

    if config.transfer_files:
        transfer_loc = config.transfer_dir
        os.makedirs(transfer_loc, exist_ok=True)

        # madgraph dir
        os.system("ls")
        print(config.process.output_dir, transfer_loc)

        copy_tree(config.process.output_dir, transfer_loc)
        
        hydra_cfg_dir = os.path.join(transfer_loc, '.hydra')
        os.makedirs(hydra_cfg_dir, exist_ok=True)
        
        print(f"hydra_cfg_dir = {hydra_cfg_dir}")

        shutil.copy('.hydra/hydra.yaml', hydra_cfg_dir)
        shutil.copy('.hydra/config.yaml', hydra_cfg_dir)
        shutil.copy('.hydra/overrides.yaml', hydra_cfg_dir)        
        
        # log files 
        [shutil.copy(f, transfer_loc) for f in glob.glob("*log*")]

        # root files
        [shutil.copy(f, transfer_loc) for f in glob.glob(".root")]
        
        