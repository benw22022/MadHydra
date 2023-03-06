import uproot
import matplotlib.pyplot as plt
import os
import glob
from pathlib import Path


def plot_results():

    files = glob.glob("rivet_output/*")

    xsec_list = []
    vevd_list = []
    mass_list = []

    for fpath in files:
        mass, vevd = fpath.split('/')[1].split('_')[-2:]

        mass_list.append(mass)
        vevd_list.append(vevd)

        xsec = 0

        if Path(f"{fpath}/output.root").is_file():

            data = uproot.open(f'{fpath}/output.root:metadata')
            xsec = data['total_xs'].array(library='np')[0]
            xsec_list.append(xsec)

        print(f"mass = {mass} GeV \t vevd = {vevd} \t xsec = {xsec}")

if __name__ == "__main__":
    plot_results()