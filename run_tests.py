"""
Testing steering script
"""

import logger
log = logger.get_logger(__name__)

import tests

if __name__ == "__main__":

    print("Running tests")
        
    tests.test_rivet()
