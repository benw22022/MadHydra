"""
Testing steering script
"""

import logging
log = logging.getLogger(__name__)
# TODO: logging like this doesn't work without hydra 

import tests

if __name__ == "__main__":

    print("Running tests")
        
    tests.test_rivet()
