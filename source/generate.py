import logger
log = logger.get_logger(__name__)

import os
import shutil
import inspect
import glob
from omegaconf import DictConfig, OmegaConf
from distutils.dir_util import copy_tree
import traceback
import source.utils as utils
import rivet
import athgeneration

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
    If requested by config will also run rivet routine provided
    Also has functionality to perform some clean up of the MadGraph production directory
    and will also transfer the files if running on batch to given folder
    args:
        config: DictConfig - config object
    returns:
        None
    """

    log.info("Running local generation")
    log.debug(f"Process is: \n {config.process}")

    # Check if we're trying to run an official dsid (we can use Athena for this)
    if config.process.get("dsid", False):
        athgeneration.generate_dsid(config)

    else:
        # Some hacky stuff to workout where MadGraph is 
        this_filename = inspect.getframeinfo(inspect.currentframe()).filename
        this_filepath = os.path.dirname(os.path.dirname(os.path.abspath(this_filename)))
        mg_path = os.path.join(this_filepath, config.madgraph_dir)
        madgraph_exec = os.path.join(mg_path, 'bin', 'mg5_aMC')


        # write the process card for MadGraph
        write_proc_card(config)
        
        # If we're running in debug mode stop here
        if config.debug: return
        
        # Run generation using local MadGraph installation
        utils.launch_process([madgraph_exec, 'proc_card.dat'], 'MadGraph', logfile='log.generate')
        
    # Run cleanup
    if config.cleanup:
        log.info("Running cleanup")
        try:
            run_cleanup(config)
        except Exception as e:
            log.error("An error occured while trying to clean up directory")
            os.system("pwd")
            os.system("ls")
            log.error(traceback.format_exc())
            log.error(e)

    # Transfer output
    if config.transfer_files:
        log.info(f"Copying file to {config.transfer_files}")
        transfer_loc = config.transfer_dir
        os.makedirs(transfer_loc, exist_ok=True)

        copy_tree(config.process.output_dir, transfer_loc)
        
        hydra_cfg_dir = os.path.join(transfer_loc, '.hydra')
        os.makedirs(hydra_cfg_dir, exist_ok=True)

        shutil.copy('.hydra/hydra.yaml', hydra_cfg_dir)
        shutil.copy('.hydra/config.yaml', hydra_cfg_dir)
        shutil.copy('.hydra/overrides.yaml', hydra_cfg_dir)        
        
        [shutil.copy(f, transfer_loc) for f in glob.glob("*log*")]
        [shutil.copy(f, transfer_loc) for f in glob.glob("*.root")]
    
    # Run rivet routine
    if not OmegaConf.is_missing(config, config.process['rivet']):
        log.info(f"Running rivet routine {config.process.rivet}")

        # If we use AthGeneration we'll get a .EVNT.root as an output, otherwise it'll be a .hepmc.gz
        file_type = '.hepmc.gz' if not config.process.get("dsid", False) else '.EVNT.root'

        rivet.rivet_analyze_job(config, file_type=file_type)

        # Copy files
        if config.transfer_files:
            [shutil.copy(f, transfer_loc) for f in glob.glob("*log*")]
            [shutil.copy(f, transfer_loc) for f in glob.glob("*.root")]

        # Also mirror to another location for convenience if requested (yeah this is kinda hacky)
        if config.get("rivet_copy_dir", False) and config.transfer_files:
            os.makedirs(config.rivet_copy_dir, exist_ok=True)
            [shutil.copy(f, config.rivet_copy_dir) for f in glob.glob("*log*")]
            [shutil.copy(f, config.rivet_copy_dir) for f in glob.glob("*.root")]

        
        
    log.info("Done")