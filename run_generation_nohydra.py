import os
import source
import argparse
from omegaconf import OmegaConf

OmegaConf.register_new_resolver("eval", eval) # This will allow us to do arithmetic in our yaml cfgs

def run_generation_nohydra():
    """
    Run generation without invoking hydra - useful for runnig generation on htcondor batch system
    Takes one mandatory command line arguement: hydra_dir
    This is the path to the .hydra director containing the config.yaml and hydra.yaml config files
    We need this to rebuild the DictConfig object and run without needing to invoke hydra
    TODO: There is probably a nicer way to do this in `run_generation.py` with hydra's compose API but I can't see it
    """
    
    # Parse arguements
    parser = argparse.ArgumentParser()
    parser.add_argument("hydra_dir", help="location of .hydra dir containing configs", type=str)
    args = parser.parse_args()

    job_conf_fpath = os.path.join(args.hydra_dir, "config.yaml")
    hydra_conf_fpath = os.path.join(args.hydra_dir, "hydra.yaml")
    
    job_conf = OmegaConf.load(job_conf_fpath)
    hydra_conf = OmegaConf.load(hydra_conf_fpath)

    config = OmegaConf.merge(job_conf, hydra_conf)

    source.run_local_generation(config)
    
if __name__ == "__main__":
    run_generation_nohydra()