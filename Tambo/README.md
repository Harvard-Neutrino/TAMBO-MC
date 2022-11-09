# Hola y bienvenidos.
This quickstart guide is the closest thing that we have to documentation of the TAMBO code, and is very much a work in progress.
Some of the functions and utilities we will go over have docstrings in the code, but many others do not.
With that in mind, if things don't work, or there is a change you would like to see, or something is just annoying but you don't know what the fix would look like, please email [Jeff Lazar](mailto:jlazar@icecube.wisc.edu).
Or if you feel up to it, you can fork this repo, make the change, and initiate a pull request.
Alright, let's get things moving.

## 0. Prerequisites

### 0.0 Jump into Julia
The code reies mostly on the Julia language so you will need to install that.
This was developed using Julia version 1.8.2 so get that one if you can.
You can find precompiled Julia version for many systems [here](https://julialang.org/downloads/).

If you are new to Julia, it isn't so bad, but there are definitely some differences from Python.
To run the code, you won't need anything advanced, but it might be good to dip you toes into some of the basics.
When I started with Julia, I had a decent handle on Python, and I used [these resources](https://github.com/JuliaAcademy/JuliaTutorials) to fill in the gaps.
There were still some growing pains, but I felt it gave me enough information to get started.

Also, I highly recommend joining the [Julia Slack](https://julialang.org/slack/).
it is very active, and everyone there is extremely helpful.

### 0.1 Install PROPOSAL
PROPOSAL is the library that we use to propoagate charged leptons.
It is written in c++, but we rely on the Python interface.
This can sometimes be installed with `pip` by running `pip install proposal==7.4.2`.
You can also install it from source using their [GitHub](https://github.com/tudo-astroparticlephysics/PROPOSAL/tree/7.4.2).
This can be a bit of battle so fortify your spirit.
If you find yourself banging your head against the wall, get in touch because we may have some tricks up our sleeve.

## 1. Conventions

### 1.0. Units

Any variable that is saved internally will *always* be in natural units with `eV=1`.
`PROPOSAL` and `CORISKA`, which we use for charged lepton propagation and air-shower simulation respectively, have different unit conventions.
These conventions are summarized in the following table.

|           | `TAMBO`    | `PROPOSAL` | `CORSIKA` |
| :---      |  :----:    |   :----:   | :---:     |
| length    | eV $^{-1}$ | cm         | cm [^1]   |
| energy    |            |            |           | 
| Footnotes | [^1] Caveat emptor | | |

If you are accessingor setting a variable that does not exist solely to go to an external dependency, *please* remember and respect this convention.

## 2. Fun

### 2.0 Download and activate the code

Use the command line to navigate to where you would like to store the repository and clone this package with:
```
git clone https://github.com/Harvard-Neutrino/TAMBO-MC.git
```
Now launch up a Julia REPL session and let's run some code.
You can activate the TAMBO code by running
```julia-repl
julia> using Pkg; Pkg.activate("TAMBO-MC/Tambo/")
  Activating project at "~/research/TAMBO-MC/Tambo"
```
And now let's import the package like
```julia-repl
julia> using Tambo
[Info: Precompiling Tambo [d9a96183-4919-46da-8188-64ea4e10e0ed]
```
Hopefully that worked...
Assuming it did, we can actually start using the code.

## 3. The `Simulator` `struct`

This is a backbone of the code, and we can make it in the following way
```julia-repl
julia> simulator = Simulator()
General configuration
_____________________
n: 10
seed: 0

Geometry configuration
______________________
geo_spline_path: /Users/jlazar/research/TAMBO-MC/resources/tambo_spline.jld2
tambo_coordinates: (-15.665°, -72.155°)


Injection configuration
_______________________
ν_pdg: 16
γ: 1.0
emin: 1.0e6 GeV
emax: 1.0e9 GeV
θmin: 0.0°
θmax: 180.0°
ϕmin: 0.0°
ϕmax: 360.0°
r_injection: 900.0 m
l_endcap: 1000.0 m
diff_xs_path: /Users/jlazar/research/TAMBO-MC/resources/cross_sections/tables/csms_differential_cdfs.h5

PROPOSAL configuration
______________________
ecut: Inf GeV
vcut: 0.01
do_interpolate: true
do_continuous: true
tablespath: /Users/jlazar/research/TAMBO-MC/resources/proposal_tables
```
We can access and set any of these fields if we want in the following way.
```julia-repl
julia> simulator.n
10

julia> simulator.n = 1_000;

julia> simulator.n
1000

julia> simulator.diff_xs_path
"/Users/jlazar/research/TAMBO-MC/resources/cross_sections/tables/csms_differential_cdfs.h5"
```
Hopefully, this last output will be different for you since you are not using my computer.
If it's not that's not great and you should change it manually, and let me know.
Also note that since we use natural units internally, some of the fields will look different when displayed in this manner.
They can be converted to your favorite unit using `units`
```julia-repl
julia> simulator.emax
1.0e18

julia> simulator.emax / units.GeV
1.0e9

julia> simulator.emax / units.MeV
1.0e12
```

Also, we can change the arguments of `Simulator` when we construct it in the following manner
```julia-repl
julia> simulator = Simulator(emin=1e5 * units.GeV)
General configuration
_____________________
n: 10
seed: 0

Geometry configuration
______________________
geo_spline_path: /Users/jlazar/research/TAMBO-MC/resources/tambo_spline.jld2
tambo_coordinates: (-15.665°, -72.155°)


Injection configuration
_______________________
ν_pdg: 16
γ: 1.0
emin: 100000.0 GeV
emax: 1.0e9 GeV
θmin: 0.0°
θmax: 180.0°
ϕmin: 0.0°
ϕmax: 360.0°
r_injection: 900.0 m
l_endcap: 1000.0 m
diff_xs_path: /Users/jlazar/research/TAMBO-MC/resources/cross_sections/tables/csms_differential_cdfs.h5

PROPOSAL configuration
______________________
ecut: Inf GeV
vcut: 0.01
do_interpolate: true
do_continuous: true
tablespath: /Users/jlazar/research/TAMBO-MC/resources/proposal_tables
```
Before we can get to simulating, I should point out there are two fields that are not listed here, `simulator.injected_events` and `simulator.proposal_events`.
Both of these are empty `Vector`s currently but that will soon.
Now, to inject some particles and the propagate the charged final products !
A warning: this will take 20 or so minutes that first time you do it since PROPOSAL will need to make tables, so I would take this time to go for a walk, grab a warm beverage of your choosing, or chitchat with a coworker.
This same 20 minutes will be needed any time you change `simulator.ecut` or `simulator.vcut` to a new value.
Once tables are made, it should take a handful of seconds to simulate 1,000 events.
```julia-repl
julia> simulator(track_progress=false);

julia> length(simulator.injected_events)
1000

julia> simulator.injected_events[1]
initial_state:
pdg_mc: 16,
energy (GeV): 4.195636776889654e6,
position (m): [129.63665506830162, -557.8187789378879, -893.3008352589561],
direction: θ (degrees): 43.59087230458852°
ϕ (degrees): 30.949511709031906°
proj: [0.5913331252113109, 0.3545999622485903, 0.7242817143909699]


final_state:
pdg_mc: 15,
energy (GeV): 4.1692218834799514e6,
position (m): [129.63665506830162, -557.8187789378879, -893.3008352589561],
direction: θ (degrees): 43.59087230458852°
ϕ (degrees): 30.949511709031906°
proj: [0.5913331252113109, 0.3545999622485903, 0.7242817143909699]
```
You can serialize this to a jld file and load it up again in the following manner
```julia-repl
julia> fname = tempname();

julia> save_simulation(simulator, fname)

julia> using JLD2; loadfile = jldopen(fname)
JLDFile /private/var/folders/4q/ncd7kk_j2t9_syxd3gmst9mc0000gn/T/jl_2hIrM9irm9 (read-only)
 ├─🔢 injected_events
 ├─🔢 proposal_events
 └─🔢 config

 julia> loadfile["config"]
Dict{Symbol, Any} with 20 entries:
  :diff_xs_path      => "/Users/jlazar/research/TAMBO-MC/resources/cross_sections/tables/csms_differential_cdfs.h5"
  :tablespath        => "/Users/jlazar/research/TAMBO-MC/resources/proposal_tables"
  :θmin              => 0.0
  :n                 => 1000
  :ϕmax              => 6.28319
  :do_interpolate    => true
  :ν_pdg             => 16
  :r_injection       => 4.56096e9
  :emax              => 1.0e18
  :do_continuous     => true
  :ϕmin              => 0.0
  :vcut              => 0.01
  :geo_spline_path   => "/Users/jlazar/research/TAMBO-MC/resources/tambo_spline.jld2"
  :tambo_coordinates => (-15.665°, -72.155°)…
  :γ                 => 1.0
  :θmax              => 3.14159
  :emin              => 1.0e14
  :l_endcap          => 5.06773e9
  :ecut              => Inf
  :seed              => 0
```
You can also construct a new simulator object from the parameters of the saved file with
```julia-repl
julia> new_simulator = Simulator(fname)
General configuration
_____________________
n: 1000
seed: 0

Geometry configuration
______________________
geo_spline_path: /Users/jlazar/research/TAMBO-MC/resources/tambo_spline.jld2
tambo_coordinates: (-15.665°, -72.155°)


Injection configuration
_______________________
ν_pdg: 16
γ: 1.0
emin: 100000.0 GeV
emax: 1.0e9 GeV
θmin: 0.0°
θmax: 180.0°
ϕmin: 0.0°
ϕmax: 360.0°
r_injection: 900.0 m
l_endcap: 1000.0 m
diff_xs_path: /Users/jlazar/research/TAMBO-MC/resources/cross_sections/tables/csms_differential_cdfs.h5

PROPOSAL configuration
______________________
ecut: Inf GeV
vcut: 0.01
do_interpolate: true
do_continuous: true
tablespath: /Users/jlazar/research/TAMBO-MC/resources/proposal_tables
```

## 4. `Simulator` under the hood

WIP

## CORSIKA interface

We are using `CORSIKA` to model the $\tau^{\pm}$ air shower.
This is a `FORTRAN 77` script that needs configuration cards to run and thus cannot be easily interfaced with our Julia code.
As a design rule, we do not store any variables internally that are not in our unit system or in the coordinate system established by the `Geometry` object.
If you are accessing or setting a variable, *please* keep this in mind.

The CORSIKA coordinate system aligns its $x-$ and $y-$ axes with the magnetic north and west respectively.
TAMBO aligns these axes with east and north.
Furthermore, the TAMBO coordinate system is always centered on the TAMBO detector, whereas CORSIKA centers it's coordinate system on the primary particles incident particle's initial horizontal position in the $x-$ and $y$- directions, and at sea-level in the $z$-direction.
CORSIKA has base units


As of this writing, we are approximating the mountain face as a perfect plane, since `CORSIKA`'s geometry handling is a bit limited.
Thus, much of the convenience functions use the internal `Plane` struct.


### To be continued...