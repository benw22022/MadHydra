"""
Plotting utilities
_____________________________________________
Utility functions for plotting results
"""

import os
import glob
import yaml
import uproot
import operator
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from functools import cache
from dataclasses import dataclass


@dataclass
class Result:
    vevd: float
    mdpp: float
    quantity: float
    file: str

    def __lt__(self, other):
        return self.vevd < other.vevd

@cache
def get_results(samples_dir, quantity):
    
    override_cgfs_paths = glob.glob(f"{samples_dir}/*/.hydra/overrides.yaml")
    
    results = []

    for fpath in override_cgfs_paths:
        with open(fpath) as f:
            cfg_list = yaml.safe_load(f)
            vevd = -999
            mdpp = -999
            for item in cfg_list:
                if 'vevD' in item:
                    vevd = 10**float(item.split("=")[-1])
                if 'MDPP' in item:
                    mdpp = float(item.split("=")[-1])

            root_file = f"{os.path.dirname(os.path.dirname(fpath))}/output.root"

            try:
                data = uproot.open(f"{root_file}:metadata")
                value = data[quantity].array(library='np')[0]
            except Exception:
                value=None

            results.append(Result(vevd, mdpp, value, root_file))

    results.sort()
    return results


def plot_results(results, name, log):

    vevds = set([r.vevd for r in sorted(results, key=operator.attrgetter('vevd'))])
    mdpps = set([r.mdpp for r in sorted(results, key=operator.attrgetter('mdpp'))])

    hmap = []

    for v in vevds:
        results_msorted = sorted(results, key=operator.attrgetter('mdpp'))
        hmap.append([r.quantity for r in results_msorted if r.vevd == v])
        
    hmap = np.array(hmap, dtype=float)

    fig, ax = plt.subplots()
    if log:
        sns.heatmap(hmap, annot=True, cmap="Oranges", xticklabels=mdpps, yticklabels=vevds, ax=ax,
                    fmt=".2", norm=LogNorm(), annot_kws={"size": 15 / np.sqrt(len(hmap))},
                    cbar_kws={'label': name});
    if not log:
        sns.heatmap(hmap, annot=True, cmap="Oranges", xticklabels=mdpps, yticklabels=vevds, ax=ax,
                    fmt=".2", annot_kws={"size": 15 / np.sqrt(len(hmap))},
                    cbar_kws={'label': name});

    ax.set_xlabel(r"$M_{\Delta^{\pm\pm}}$ [GeV]");
    ax.set_ylabel(r"vevD [GeV]");