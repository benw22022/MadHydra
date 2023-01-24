"""
Test rivet routine running
"""

import logging
log = logging.getLogger(__name__)

import os
import rivet
from omegaconf import OmegaConf

def test_rivet() -> None: 
    
    TESTLOC = 'generate_output/2022-11-30_16-54-18_typeII_single_prod_WWjj_502_1e-07'
    
    test_conf = OmegaConf.load(os.path.join(TESTLOC, '.hydra', 'config.yaml'))
    hydra_test_conf = OmegaConf.load(os.path.join(TESTLOC, '.hydra', 'hydra.yaml'))
    test_conf.merge_with(hydra_test_conf)
    
    print("Running rivet test")
    
    rivet.run_rivet_routine(test_conf)