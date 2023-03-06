import logging
log = logging.getLogger(__name__)

import os
import htcondor
import classad
from omegaconf import DictConfig
from hydra.utils import get_original_cwd, to_absolute_path


def submit_job(config: DictConfig):

    cmd = f"#!/bin/bash \n\
eval \"$(conda shell.bash hook)\" \n\
conda activate {config.conda_env} \n\
python3 {get_original_cwd()}/run_generation_nohydra.py {os.getcwd()}/.hydra \n\
"

    with open("submit.sh", 'w') as file:
        file.write(cmd)

    config.condor.executable = 'submit.sh'

    if config.debug:
        log.info("Debug mode set - won't actually submit job")
        return

    hostname_job = htcondor.Submit(config.condor)

    schedd = htcondor.Schedd()                   
    submit_result = schedd.submit(hostname_job)  
    print(submit_result.cluster())               

