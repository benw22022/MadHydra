import logging
log = logging.getLogger(__name__)

import subprocess
from subprocess import Popen, PIPE, STDOUT
from typing import List
import re

def escape_ansi(line):
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

def launch_process(cmd_list: List[str], proc_name: str='', shell: bool=False) -> int:
    """
    Launches a subprocess, just a little function to avoid some of the boilerplate
    args:
        cmd_list: List[str] - a list of commands/options to execute. e.g "rm -r aDir" would be ['rm', '-r', 'aDir']
        proc_name: str (default='')- a string to optionally label the process with when logging
    returns:
        int: exit code (0 means success)
    """
    
    if len(cmd_list) == 0:
        return 0
    
    process = Popen(cmd_list, stdout=PIPE, stderr=STDOUT, shell=shell)
    with process.stdout:
        log_subprocess_output(process.stdout, proc_name)
        
    retcode = process.wait()
        
    return retcode