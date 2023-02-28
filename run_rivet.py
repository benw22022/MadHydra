import os
import glob
import argparse
from omegaconf import OmegaConf
from rivet import rivet_analyze_job
from multiprocessing import Pool

def run_rivet_in_new_dir(cfg, routine, run_folder='rivet_output'):
    
    cwd = os.getcwd()
    run_dir = os.path.join(run_folder, cfg.process.output_dir)
    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)

    rivet_analyze_job(cfg, routine=routine)

    os.chdir(cwd)

def main() -> None:

    # Parse arguements
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", help="directory to look through for simulation jobs", type=str)
    parser.add_argument("--routine", help="if given, overides rivet routine provided in original config", type=str, default='')
    parser.add_argument("--mp", help="Use multiprocessing", type=bool, default=False)
    args = parser.parse_args()

    print(args.routine)

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
    if args.mp:
        pool = Pool()
        for cfg in configs:
            pool.apply_async(run_rivet_in_new_dir, args = (cfg, args.routine))
            pool.close()
            pool.join()
    
    else:   
        for cfg in configs:
            run_rivet_in_new_dir(cfg, args.routine)


if __name__ == "__main__":
    main()
