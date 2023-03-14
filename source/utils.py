import logger
log = logger.get_logger(__name__)

import os
import re
from pathlib import Path
from typing import List
from subprocess import Popen, PIPE, STDOUT


def escape_ansi(line: str):
    """
    An attempt to clense stdout of ansi escape charaters for logfile writing
    TODO: This doesn't really work - improvement welcome!
    Args:
        line (str): line from stdout
    Returns:
        str: A new line with (hopefully) the ansi escape characters removed
    """
    try:
        ansi_escape =re.compile(r'(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]')
        return ansi_escape.sub('', line)
    except TypeError:
        return f"Unable to parse {line}"
        

def log_subprocess_output(pipe, proc_name: str='') -> None:
    """
    Logs output from subprocess to log file
    args:
        pipe: stdout buffer to read from
        proc_name: str (default='')- a string to optionally label the process with
    returns:
        None
    """
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        cleaned_line = escape_ansi(str(line))
        cleaned_line = str(cleaned_line).replace(r"\n", "").strip().replace("b\"", "").rstrip("\'").lstrip("\'")
        log.info(f"{proc_name}: %s", cleaned_line)

def launch_process(cmd_list: List[str], proc_name: str='', shell: bool=False, env_vars=None, logfile: str=None) -> int:
    """
    Launches a subprocess, just a little function to avoid some of the boilerplate
    args:
        cmd_list: List[str] - a list of commands/options to execute. e.g "rm -r dirname" would be ['rm', '-r', 'dirname']
        proc_name: str (default='') - a string to optionally label the process with when logging
        shell: bool (default=False) - If True then use Popen(..., shell=True)
    returns:
        int: exit code (0 means success)
    """
    
    my_env = os.environ.copy()

    if env_vars is not None:
        my_env = {**my_env, **env_vars}

    if len(cmd_list) == 0:
        log.warnng(f"tried to launch subprocess ({proc_name}) with no arguements - returning `0` instead")
        return 0
    
    stdout=PIPE
    # if logfile: # TODO: workout how to redirect to log file
    #     stdout = open(logfile, 'w')

    process = Popen(cmd_list, stdout=stdout, stderr=STDOUT, shell=shell, env=my_env)

    with process.stdout:
        log_subprocess_output(process.stdout, proc_name)
        
    retcode = process.wait()
    
    # if logfile:
    #     stdout.close()

    return retcode


def delete_files_with_extn(directory: str, extn: str) -> List[str]:
    """
    Recursivly searches though folders to delete files with matching extension
    args:
        directory: str - path to search through
        extn: str - file extension to delete
    returns:
        deleted_files: List[str] - list of files that had been deleted
    """

    deleted_files = []

    for path in Path(directory).rglob(f'*{extn}'):
        os.remove(path)
        log.debug(f"Deleted file: {path}")
        deleted_files.append(path)

    return deleted_files


def get_files_with_extn(directory: str, extn:str) -> List[str]:
    """
    Get all files in directory with exension
    args:
        directory: str - directory to search
        extn: str - file extension to search for
    returns:
        List[str] - list of file paths of files ending in `extn`
    """
    return [os.path.join(r, fn)
        for r, ds, fs in os.walk(directory) 
        for fn in fs if fn.endswith(extn)]