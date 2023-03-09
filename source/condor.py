import logger
log = logger.get_logger(__name__)

import os
import htcondor
from omegaconf import DictConfig
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


    cmd = f"#!/bin/bash \n\
# eval \"$(conda shell.bash hook)\" # You can use this on lxplus\n\
# conda activate {config.conda_env} \n\
source  {config.conda_dir}/bin/activate {config.conda_env} \n\
python3 {get_original_cwd()}/run_generation_nohydra.py {os.getcwd()}/.hydra \n\
"

    with open("submit.sh", 'w') as file:
        file.write(cmd)

    config.condor.executable = "submit.sh"

    if config.debug:
        log.info("Debug mode set - won't actually submit job")
        return

    with open("htc.submit", 'w') as file:
        for key, value in config.condor.items():
            file.write(f"{key} = {value}\n")
        file.write("queue")

    # os.system("condor_submit htc.submit")
    hostname_job = htcondor.Submit(config.condor)

    schedd = htcondor.Schedd()                   
    submit_result = schedd.submit(hostname_job)   # TODO Add back in batch-name to jobs
    log.info(f"submitted job {submit_result.cluster()}")               

