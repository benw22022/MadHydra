"""
Test rivet routine running
"""

import logging
log = logging.getLogger(__name__)
# TODO: setup logging without Hydra

import os
import rivet
from omegaconf import OmegaConf

def test_rivet() -> None: 
    
    # Test data
    TESTLOC = 'tests/data/2023-01-24_13-39-16_typeII_single_prod_WWjj_100_1.0'
    
    # Join configs
    test_conf = OmegaConf.load(os.path.join(TESTLOC, '.hydra', 'config.yaml'))
    hydra_test_conf = OmegaConf.load(os.path.join(TESTLOC, '.hydra', 'hydra.yaml'))
    test_conf.merge_with(hydra_test_conf)
    
    # Move to run dir
    os.chdir(test_conf.hydra.runtime.output_dir)

    print("Running rivet test")
    
    rivet.run_rivet_routine(test_conf, keep_gz_files=True)