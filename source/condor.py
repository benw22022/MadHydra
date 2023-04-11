import logger
log = logger.get_logger(__name__)

import os
import htcondor
from omegaconf import DictConfig, OmegaConf
from hydra.utils import get_original_cwd


def submit_job(config: DictConfig) -> None:
    """
    Uses the htcondor python bindings to submit a job to condor batch system
    Note: Had issues with using the htcondor python API on lxplus - so you may need to use os.system instead
    args:
        config: DictConf - config object to generate
    returns:
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

    hostname_job['arguments'] = arguments   
    hostname_job['job_name'] = arguments

    if config.debug:
        log.info(f"Job arguments are {arguments}")
        log.info("Debug mode set - won't actually submit job")
        return

    # now submit job
    schedd = htcondor.Schedd()                   
    submit_result = schedd.submit(hostname_job)   # TODO Add back in batch-name to jobs
    log.info(f"submitted job {submit_result.cluster()}: {arguments}")               

