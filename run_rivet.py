import os
import glob
import argparse
from omegaconf import OmegaConf
from rivet import rivet_analyze_job


def main() -> None:

    # Parse arguements
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", help="directory to look through for simulation jobs", type=str)
    parser.add_argument("routine", help="if given, overides rivet routine provided in original config", type=str, default='')
    args = parser.parse_args()

    # Get .hydra dirs, check if we're running over a single simulation job. Otherwise look onelayer deeper
    hydra_dirs = glob.glob(os.path.join(args.run_dir, ".hydra"))
    if len(hydra_dirs) == 0:
        hydra_dirs = glob.glob(os.path.join(args.run_dir, "*", ".hydra"))

    # Load and merge job & hydra configs
    configs = []
    for conf_dir in hydra_dirs:
        job_conf_fpath = os.path.join(conf_dir, "config.yaml")
        hydra_conf_fpath = os.path.join(conf_dir, "hydra.yaml")
        
        job_conf = OmegaConf.load(job_conf_fpath)
        hydra_conf = OmegaConf.load(hydra_conf_fpath)

        configs.append(OmegaConf.merge(job_conf, hydra_conf))

    # Loop through configs and run routines
    for cfg in configs:
        rivet_analyze_job(cfg, routine=args.routine)

if __name__ == "__main__":
    main()
