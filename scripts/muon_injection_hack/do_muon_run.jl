project_dir = (@__DIR__) * "/../../"
using Pkg
Pkg.activate(project_dir)
using Tambo
using StaticArrays

spline_path = project_dir * "resources/perfect_valley_spline.jld2"
tambo_coordinates = Coord( deg2rad.([0.0, -0.015738128103577786])... )
normal_vec = Tambo.Direction( [0.5655283556, 0, 0.8247288518]... )

geo = Geometry(spline_path, tambo_coordinates, normal_vec)

pdg = 13 # muon
energy = 1.0e3 * units.TeV
zenith = 90.0 * units.degree
azimuth = 180.0 * units.degree
inject_pos = SVector{3}([3.8191921355258502, 0.0004210065329517819, 0.00023716158800735456]) * units.km
intercept_pos = SVector{3}([-0.00034386310581669354, 0.00042026995597829026, 0.00023579184402179322]) * units.km
plane = SVector{3}([0.565528355605942, 0.0, 0.8247288518086653])
obs_z = 2.851251161459321 * units.km
thinning = 1.0e-6
ecuts = SVector{4}([0.001, 0.001, 0.05, 0.05] * units.GeV) # em, photon, mu, hadron cuts
corsika_path = "/home/lordkelvin/packages/corsika/corsika-install/bin/c8_air_shower"
corsika_FLUPRO = "/home/lordkelvin/packages/fluka"
corsika_FLUFOR = "gfortransbatch"
outdir = "./"
proposal_index = 1
decay_index = 1
seed = 0
parallelize_corsika = false

Tambo.corsika_run(
    pdg,
    energy,
    zenith,
    azimuth,
    inject_pos,
    intercept_pos,
    plane,
    obs_z,
    thinning,
    ecuts,
    corsika_path,
    corsika_FLUPRO,
    corsika_FLUFOR,
    outdir,
    proposal_index,
    decay_index,
    seed,
    parallelize_corsika=parallelize_corsika
    )