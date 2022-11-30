MadHydra
========================

This repository acts as a framework for running [`MadGraph`](https://launchpad.net/mg5amcnlo) through the use of yaml config files using the [`hydra`](https://hydra.cc/) library. Though `hydra` we can quickly and easily generate new directory structures, organise processes and parameters in config files and easily execute multirun jobs. Hopefully, this should make your studies easier to organise and run.
Although primarily used for configuring machine learning projects, `hydra` can also be used to configure MadGraph generation jobs.
Note: This documentation does not intend to be a tutorial for both `hydra` and `MadGraph`

Installation
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

Once you have cloned the repository and setup your python enviroment, you'll need to setup `MadGraph`. Go to the `MadGraph` [website](https://launchpad.net/mg5amcnlo) and download the latest version. This will give you a compressed tar file. Copy the file into this repository and untar it with:
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

Generating events
--------------------------

To run event generation in `MadHydra` you'll need to define your own process config file. This is a `yaml` file placed in `config/processes` and look something like this template:
```
---
# Name of model
model: MyBSMModel 

# Name of process generated
process_name: NameOfMyProcess

# Name of directory that MadGraph will create where it will run the generation
output_dir: ${process.process_name}_${parameters.param1}_${parameters.param2}

# Ordered list of commands to give to MadGraph
gen_cmd: ['import model ${process.model}',
          'generate p p > x y z',
          'output ${process.output_dir}',
          'launch ',
          'shower=pythia8',
          'set param1: ${parameters.param1}',
          'set param2: ${parameters.param2}'
          ]
```

Notice how you can set placeholders for parts of strings using the `${<myConfigItem>}` syntax. 
For our model parameters we set these in a differnt config file in `config/parameters.yaml`, which contains the `parameters` field:
```
--- 
parameters:
    param1: 100 
    param2: 0.001
```

To run the event generation all you need to do is tell `hydra` which config you want to use for process

```
python3 run_generation.py +process=<my_process>
```
Where `<my_process>` is simply the name of your process config file without the `.yaml`. e.g if you have a config called `SM_ssWW_EWK.yaml` you would generate it with the command
```
python3 run_generation.py +process=SM_ssWW_EWK
```

When you run this command, a new directory will be created here `generate_output/<year><month><day><hour><min><sec>_SM_ssWW_EWK` for `MadGraph` to run in.

You'll also notice that inside this directory we will also get a `.hydra` directory, inside you'll find the configuration that were used to run this job. Never lose track of what settings you used to generate your events ever again!

Running on htc batch system
------------------------------------
If you want to run on the htc condor batch system you can, all you have to do is set `batch=True` e.g.
```
python3 run_generation.py +process=SM_ssWW_EWK batch=True
```