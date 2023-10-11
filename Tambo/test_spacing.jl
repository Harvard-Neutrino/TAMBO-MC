using Revise
using Pkg
Pkg.activate(".")
using Tambo
using PyCall
using JLD2
using CSV
using StaticArrays

np = pyimport("numpy")
ak = pyimport("awkward")

pavel_sim = jldopen("/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.jld2")
config = SimulationConfig(; pavel_sim["config"]...)
injector = Tambo.Injector(config)
geo = Tambo.Geometry(config)
plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)

# This is the parameter we want to vary
ℓ = parse(Float64, ARGS[1]) * units.m

# This is the spacing between detectors
Δs = parse(Float64, ARGS[2]) * units.m # meters

detection_modules = Tambo.make_trianglearray(-900units.m, 2000units.m, -ℓ/2, ℓ/2, Δs, ϕ=whitepaper_normal_vec.ϕ)

# Remove detectors that are out of range
zmin = -1100units.m
zmax = 1100units.m
mask = zmin .< Tambo.plane_z.(getfield.(detection_modules, :x), getfield.(detection_modules, :y), Ref(plane)) .< zmax
detection_modules = detection_modules[mask]

outdir = ARGS[3]

neutrino_events = np.load("/n/home12/jlazar/CORSIKA_runs.npy")
basedir = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/whitepaper_sim/"
pavel_magic_csv = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.csv"
savefile = "$(outdir)/test_events_$(ℓ/units.m)_$(Δs / units.m).npy"

csv = CSV.File(pavel_magic_csv)
idxs = getindex.(csv, 16)
subidxs = getindex.(csv, 17)

function load(files, csv_x, csv_y, csy_z)
    ycorsika = SVector{3}([0.89192975455881607, 0.18563051261662877, -0.41231374670066206])
    xcorsika = SVector{3}([0, -0.91184756344828699, -0.41052895273466672])
    zcorsika = whitepaper_normal_vec.proj
    xyzcorsika = [
      xcorsika.x xcorsika.y xcorsika.z;
      ycorsika.x ycorsika.y ycorsika.z;
      zcorsika.x zcorsika.y zcorsika.z;
    ]
    particle_xs = []
    particle_ys = []
    particle_zs = []
    particle_ts = []
    weights = []
    for file in files
        x = nothing
        try
            x = ak.from_parquet(file)
        catch
            continue
        end
        xs = x["x"].to_numpy() .* units.m
        ys = x["y"].to_numpy() .* units.m
        zs = x["z"].to_numpy() .* units.m
        ws = x["weight"].to_numpy()
        new_poss = []
        for (x, y, z) in zip(xs, ys, zs)
            push!(new_poss, transpose([x,y,z]' * xyzcorsika) .+ [csv_x, csv_y, csy_z])
        end
        ts = x["time"].to_numpy() .* units.second
        particle_xs = cat(particle_xs, getindex.(new_poss, 1), dims=1)
        particle_ys = cat(particle_ys, getindex.(new_poss, 2), dims=1)
        particle_zs = cat(particle_zs, getindex.(new_poss, 3), dims=1)
        particle_ts = cat(particle_ts, ts, dims=1)
        weights = cat(weights, ws, dims=1)
    end
    return particle_xs, particle_ys, particle_zs, particle_ts, weights
end

thresh = 1units.m # distance from detection module center to consider a hit
triggered_evts = []
for (kdx, neutrino_event) in enumerate(neutrino_events)
  # Find all showers for this events
  files = Tambo.find_extant_files(neutrino_event, basedir)
  # Make a vactor that we will fill with hits
  desc = split(first(files), "/")[end-2]
  idx = parse(Int, split(desc, "_")[end-1])
  row = first(csv[idx.==idxs])
  csv_x, csv_y, csv_z = row.inter_x, row.inter_y, row.inter_z

  xs, ys, _, _, weights = load(files, csv_x, csv_y, csv_z)

  # Find hits that fall within a certain distance of any detection unit
  hits = Tambo.find_near_hits(
    xs,
    ys,
    weights,
    detection_modules;
    thresh=thresh
  )
  # A dictionary whose keys are detection modules and whose values are
  # integer number of particles. This does Poisson sampling for events
  # that have thinning
  d = Tambo.make_hit_dict(hits)
  # Count the number of particles from triggered modules
  n = 0
  for (_, v) in d
    # a module is not triggered if it has less than three hits
    if v < 3
      continue
    end
    n += v
  end
  # Trigger criterion is 30 particles from tirggered detection modules
  if n >= 20
    push!(triggered_evts, neutrino_event)
  end
  if kdx % 100==0
    println((kdx,length(triggered_evts) / kdx))
  end
end

events = pavel_sim["injected_events"][triggered_evts]
# Northern tracks stuff
γ = 2.37
norm = 1.44e-18 / units.GeV / units.cm^2 / units.second * 1e14^γ

# HESE stuff. Factor of three is becasue they report all flavor flux
# γ = 2.87
# norm = 6.37e-18 / 3 / units.GeV / units.cm^2 / units.second * 1e14^γ

pl = Tambo.PowerLaw(γ, 100units.GeV, 1e9units.GeV, norm)
fluxes = pl.(getfield.(getfield.(events, :initial_state), :energy))

weights = Tambo.oneweight.(events, Ref(injector.xs), Ref(injector.xs), Ref(injector.powerlaw), Ref(injector.anglesampler), Ref(injector.injectionvolume), Ref(geo))
rates = fluxes .* weights ./ maximum(neutrino_events)
println(sum(rates) * 10^7.5 * units.second)

output = [triggered_evts, rates]

np.save(savefile, output)
