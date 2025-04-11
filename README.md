# Hola y bienvenidos.
Welcome to TAMBOSim! This quickstart guide will show you how to install TAMBOSim and run your first TAMBO simulation. Some of the functions and utilities we will go over have docstrings in the code, but many others do not. With that in mind, if things don't work, or there is a change you would like to see, or something is just annoying but you don't know what the fix would look like, please feel free to create an [issue in the GitHub repo](https://github.com/Harvard-Neutrino/TAMBO-MC/issues). You should also feel free to email [Jeff Lazar](mailto:jlazar@icecube.wisc.edu) or [Will Thompson](mailto:will_thompson@g.harvard.edu). Lastly, if you feel up to it, you can fork this repo, make the change, and initiate a pull request yourself.
Alright, let's get things moving!

## [0] Prerequisite: Jump into Julia
The code reies mostly on the Julia language so you will need to install that. This was developed using Julia version 1.11.1. Get that one if you can, but more recent versions should work fine. You can find precompiled Julia version for many systems [here](https://julialang.org/downloads/). While this will by default install the most recent version of Julia, installing Julia in this way will also install [Juliaup](https://github.com/JuliaLang/juliaup), a Julia version manager. Version 1.11.1 can then be installed with the command `juliaup add 1.11.1`; you can make this version the default that runs when calling `julia` with the command `juliaup default 1.11.1`

If you are new to Julia, it isn't so bad, but there are definitely some differences from Python. To run the code, you won't need anything advanced, but it might be good to dip you toes into some of the basics. When I started with Julia, I had a decent handle on Python, and I used [these resources](https://github.com/JuliaAcademy/JuliaTutorials) to fill in the gaps. There were still some growing pains, but I felt it gave me enough information to get started.
Also, I highly recommend joining the [Julia Slack](https://julialang.org/slack/). It is very active, and everyone there is extremely helpful.

## [1] Installing Dependencies
The physics of TAMBOSim relies primarily on three external software packages: PROPOSAL, TauRunner, and CORSIKA. PROPOSAL and TauRunner can both be run using Python. Of course we are using Julia, not Python, so we will interface with these packages using the Julia package PyCall. CORSIKA, in contrast, is written in C++. TAMBOSim interfaces directly with the CORSIKA executable, so it must be compiled directly, which will be covered later in this section.

We have found that the most straightforward way to get the Python dependencies to play nice in TAMBOSim is by setting up a clean Python virtual environment, or venv. Let’s do that first.

### [1.1] Python venv
We’ll assume that you have a relatively recent version of Python installed. We have found that version 3.12.4 works well. Create a fresh venv by running `python -m venv /path/to/tambo_venv`. Go ahead and activate this venv using `source /path/to/tambo_venv/bin/activate`.

### [1.2] PROPOSAL
PROPOSAL is the library that we use to propoagate charged leptons. It is written in C++, but we rely on the Python interface. We’ve found installing PROPOSAL can sometimes be a bit of a battle, so fortify your spirit. If you find yourself banging your head against the wall, get in touch because we may have some tricks up our sleeve.

PROPOSAL is easiest to install using pip. Running `pip install proposal==7.4.2` to install the version used by TAMBOSim. Alternatively, you can install it from source using their [GitHub](https://github.com/tudo-astroparticlephysics/PROPOSAL/tree/7.4.2).

### [1.3] TauRunner
TauRunner is used to propagate high-energy tau neturinos thought the Earth, taking into account the effects of [tau regeneration](https://doi.org/10.48550/arXiv.hep-ph/9804354).  To install it, you will need to directly clone the TauRunner repo. After cloning the repo, install TauRunner by running `pip install /path/to/TauRunner`.

You can test that the TauRunner install was successful by running `import taurunner`. Note that the first time you actually propagate a tau PROPOSAL will first need to build its cross section tables. This can take several minutes, so don’t worry if there’s a delay.

### [1.4] PyCall
Now you can work on getting these Python packages to play nice with Julia. Go ahead and clone TAMBOSim from the repo on GitHub. We’ll assume you clone it into `/path/to/TAMBOSim/`.

Start an interactive Julia session with `julia` and activate the TAMBOSim project environment. This is done using the Pkg program with the following commands: `using Pkg, Pkg.activate(“/path/to/TAMBOSim/Tambo”)`. 

To get PyCall working, first run `ENV["PYTHON"]=/path/to/tambo_env/bin/python` (you can also get this path by running `which python` on the command line while inside your TAMBO venv). This tells PyCall which Python executable it should use. Now run `using Tambo; Pkg.build("PyCall")` to build PyCall. After this completes, try running `using PyCall; tr = pyimport("taurunner")`. If this succeeds, the installation was successful!

Note: when setting the `PYTHON` Julia environmental variable, `~` is not automatically expanded into your home directory. This means if you do not manually replace `~` with the path to your home directory, PyCall will fail to find your Python install. For I so love the users of TAMBOSim, that I sacrificed hours of my life figuring out this excentricity so that you shall not suffer as I did.

## 2. Running the Code
 ### [2.1] Configuring Your Environment
After downloading and seting up the required packages as described above, we now need to configure your environment. If you haven’t already, activate your Python venv using `source /path/to/tambo_venv/bin/activate`.

You should also define environemental variables that specify where `TAMBOSim` is installed and where your top-level data directory for all `TAMBOSim` simulations is located. In a shell, run
```
export TAMBOSIM_PATH=/path/to/TAMBOSim/Tambo
export TAMBO_DATA_PATH=/path/to/TAMBOSim/data/TAMBO_data
```
We also need to tell `TAMBOSim` where to find `CORSIKA` and files needed by `CORSIKA`. This information is handled within a simulation configuration file. Open up the config file at `$TAMBOSIM_PATH/resources/configuration_examples/test_TAMBOSim.toml` and edit
* `corsika_path` to point to your `CORSIKA8` executable
* `FLUPRO` to point to the version of `FLUKA` for `CORSIKA` to use
* `FLUFOR` to point to the version of `FORTRAN` used by `FLUKA`

### [2.2] Precompile `TAMBOSim`
Now launch up a Julia REPL session and activate the TAMBO code by running
```julia-repl
julia> using Pkg; Pkg.activate(ENV["TAMBOSIM_PATH"])
  Activating project at "/path/to/TAMBO-MC/Tambo"
```
And now let's import the package like. If this is your first time using `TAMBOSim` or some updates have been made since its last import, Julia will precompile the package on import. To import `TAMBOSim` run
```julia-repl
julia> using Tambo
[Info: Precompiling Tambo [d9a96183-4919-46da-8188-64ea4e10e0ed]
```
Hopefully that worked...
Assuming it did, we can actually start using the code.

### [2.3] Let’s Simulate Some Events
Several examples on how to use `TAMBOSim` can be found in the `examples/` directory. For now, we will run a simulation that passes through the whole simulation chain to ensure everything is running properly.

`TAMBOSim` makes use of the `Snakemake` workflow management system. More details on `Snakemake` can be found [here](https://snakemake.github.io). By default, `TAMBOSim` will output the results of the simulation to a directory with the name of the config file in the `TAMBO_DATA_PATH` directory. For example, the name of the config file we will use here is `test_TAMBOSim.cfg`, so the output files will be written to `TAMBO_DATA_PATH/test_TAMBOSim/`. Go ahead and create and change to this directory.

To do a full run of `TAMBOSim`, execute `snakemake -p -s $TAMBOSIM_PATH/scripts/Snakefile $TAMBO_DATA_PATH/test_TAMBOSim/triggered/2300_3700_5000_150_small/whitepaper_trigger/triggered_events_00000.jld2`. This will launch the full simulation chain, where `TAMBOSim` will:
* Inject tau neutrinos into the simulation volume using `TauRunner` and compute the location of the associated tau lepton decay with `PROPOSAL`
* Simulate the tau-induced air shower using `CORSIKA`
* Identify events that trigger the detector array
If all has worked as expected, Snakemake will issue you a success message and the final triggered event file will be produced.

If you run into any problems along the way, please don’t hesitate to get in contact with us. We’d love to help you get your simulation running.

### [2.4] Simulating a Lot of Events
By default, `Snakemake` can also be used to generate a `TAMBOSim`  dataset in parallel on a cluster with minimal alterations to the above procedure. All you need to do is to open up the Snakefile at `$TAMBOSIM_PATH/scripts/Snakefile` and change the variable `slurm_partition` to a comma-separated list of the slurm partitions you’d like to run the simulation on. Then, simply run `snakemake -p --executor slurm --jobs N -s $TAMBOSIM_PATH/scripts/Snakefile $TAMBO_DATA_PATH/test_TAMBOSim/triggered/2300_3700_5000_150_small/whitepaper_trigger/triggered_events_00000.jld2`, where `N` is the maximum number of slurm jobs to run in parallel.

## Citation
Tell people how to cite us.

