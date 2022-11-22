import logging
log = logging.getLogger(__name__)

from omegaconf import DictConfig

def write_proc_card(config : DictConfig) -> None:
    
    with open('proc_card.dat', 'w') as proc_card:
        for line in config.process.gen_cmd:
            proc_card.write(f'{line}\n')