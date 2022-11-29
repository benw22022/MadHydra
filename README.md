MadHydra
========================

This repository acts as a framework for running MadGraph through the use of yaml config files using the [`hydra`](https://hydra.cc/) library. 
Although primarily used for configuring machine learning projects, `hydra` can also be used to configure MadGraph generation jobs.
Note: This documentation does not intend to be a tutorial for both `hydra` and `MadGraph`

Installation and setup
-------------------
To install `MadHydra` do:
```
git clone https://github.com/benw22022/MadHydra
```
You'll then need to setup your `python` enviroment, I recommend using [`miniconda`](https://docs.conda.io/en/latest/miniconda.html). If you haven't used `miniconda` before, it allows you to setup and use self contained python enviroments which helps you avoid dependancy conflicts between different packages.
Once you've followed the installation instructions you can then create your `python` enviroment with
```
conda create -n madhydra python=3.9
```
Then activate do
```
conda activate madhydra
```
Then to install `hydra` do 
```
pip install hydra-core --upgrade
```
Optionally, you can also install the `colorlog` plugin to make your teminal messages colourful
```
pip install hydra_colorlog --upgrade
``` 

Once you have cloned the repository and setup your python enviroment, you'll need to setup `MadGraph`. Go to the `MadGraph` (website)[https://launchpad.net/mg5amcnlo] and download the latest version. This will give you a compressed tar file. Copy the file into this repository and untar it with:
```
tar -xzf MG5_aMC_v3.X.X.tar.gz
```
Where `MG5_aMC_v3.X.X.tar.gz` is the name of the tar file you downloaded.

To complete the `MadGraph` installation you'll probably want to install `pythia8` and `lhapdf6`. To do this start the `MadGraph` command line interface (CLI) do:
```
./MG5_aMC_v3_4_1/bin/mg5_aMC
```
(You'll need your base python enviroment be 3.7 or above for this to work, so if this doesn't run actiavte the conda enviroment from earlier)

Once you have opened the `MadGraph` CLI you can then do:
```
install pythia8
install lhapdf6
```

Finally, if you have a BSM physics model which you wish to generate events with, you'll need to place it in `MG5_aMC_v4_4_1/models`.

With that done, the setup of this framework is now complete. 





