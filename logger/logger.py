import logging
import sys
APP_LOGGER_NAME = 'MadHydra'

def setup_applevel_logger(logger_name = APP_LOGGER_NAME, file_name=None): 

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s][%(name)s][%(levelname)s] - %(message)s")
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)
    if file_name:
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    return logger


def get_logger(module_name):
    """
    Gets the logger
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name) 


def setup_default_logger(file_name=None):
    """
    Setup a default logger for the application if we are not runnning a scipt launched with @hydra
    """
    setup_applevel_logger(APP_LOGGER_NAME, file_name=file_name)