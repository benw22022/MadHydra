import logger
log = logger.get_logger(__name__)

import os
import htcondor
from omegaconf import DictConfig, OmegaConf
from hydra.utils import get_original_cwd
from typing import List


def write_generation_executable(config: DictConfig, arguments: List[str]) -> None:
    """
    Write bash script to source python environment and run generation
    Args:
        config (DictConfig): config 
        arguments (List[str]): python args - i.e. /home/usr/MadHydra/run_generation.py +process=...
    Returns:
        None
    """
    with open("generate.sh", 'w') as f:
        cmd = ""
        cmd += "#!/bin/bash\n"
        cmd += f"source  {config.conda_dir}/bin/activate {config.conda_env} \n"
        cmd += f"python3 {arguments} \n"
        f.write(cmd)



def submit_job(config: DictConfig) -> None:
    """
    Uses the htcondor python bindings to submit a job to condor batch system
    Note: Had issues with using the htcondor python API on lxplus - so you may need to use os.system instead
    Args:
        config (DictConf): config object to generate
    Returns:
        None
    """

    # create our condor job object
    hostname_job = htcondor.Submit(config.condor)

    # parse arguments
    overrides = OmegaConf.load(".hydra/overrides.yaml")
    run_gen_py = os.path.join(get_original_cwd(), "run_generation.py")  # path to run_generation.py

    arguments = f"{run_gen_py}"
    for override in overrides:
        if 'batch' in override:
            continue
        arguments += f" {override}"
    arguments += "__on_batch=True"

    # hostname_job['arguments'] = arguments   
    # hostname_job['job_name'] = arguments

    write_generation_executable(config, arguments)

    if config.debug:
        log.info(f"Job arguments are {arguments}")
        log.info("Debug mode set - won't actually submit job")
        return

    log.info(hostname_job)

    # return 
    # now submit job
    schedd = htcondor.Schedd()                   
    submit_result = schedd.submit(hostname_job)   # TODO Add back in batch-name to jobs
    log.info(f"submitted job {submit_result.cluster()}: {arguments}")               


