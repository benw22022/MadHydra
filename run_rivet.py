import os
import glob
import argparse
from omegaconf import OmegaConf, DictConfig
from rivet import rivet_analyze_job
from multiprocessing import Pool

def run_rivet_in_new_dir(cfg: DictConfig, routine: str, run_folder: str='rivet_output') -> None:
    """
    Creates a new directory, cds into it, runs rivet and cds out
    args:
        cfg: DictConfig - config object
        routine: str - name of rivet routine to run (must be in rivet/routines). Does not need .cc extsn
        run_folder: str (default='rivet_output') - location to store results
    returns:
        None
    """
    cwd = os.getcwd()
    run_dir = os.path.join(run_folder, cfg.process.output_dir)
    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)

    rivet_analyze_job(cfg, routine=routine)

    os.chdir(cwd)


def main() -> None:
    """
    Standalone, non-hydra function for running rivet on MadHydra jobs post-generation
    Useful if for some reason rivet was not scheduled to run or you want to reanalyse the data with a new routine
    Script takes one compulsory arguement - 'run_dir' - this is the location of either
    a standalone MadHydra generation folder or a directory containing such folders.
    In either case there must be a .hydra directory in the 1st directory level, for single run case; or 
    2nd level in multirun case.
    E.g. Single run                     Multirun   
        - <process-name>                - <multirun-dir>        
            - .hydra                       - <process-name-1>
            - <process-name>                   - .hydra
                                               - <process-name-1>        
                                           - <process-name-2>
                                               - .hydra
                                               - <process-name-2>        
    The other arguements are:
        --routine: here you can override the existing routine in the .hydra/config.yaml file
        --mp: if true, use python multiprocessing. Will run rivet analysis jobs in parallel
    """

    # Parse arguements
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", help="directory to look through for simulation jobs", type=str)
    parser.add_argument("--routine", help="if given, overides rivet routine provided in original config", type=str, default='')
    parser.add_argument("--mp", help="Use multiprocessing", type=bool, default=False)
    parser.add_argument("--cores", help="number of cores to use for mp", type=bool, default=None)
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
        pool = Pool(args.cores)
        for cfg in configs:
            pool.apply_async(run_rivet_in_new_dir, args = (cfg, args.routine))
        pool.close()
        pool.join()

    else:   
        [run_rivet_in_new_dir(cfg, args.routine) for cfg in configs]


if __name__ == "__main__":
    main()
