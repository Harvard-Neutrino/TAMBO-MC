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

const zmin = -1100units.m
const zmax = 1100units.m
const ycorsika = SVector{3}([0.89192975455881607, 0.18563051261662877, -0.41231374670066206])
const xcorsika = SVector{3}([0, -0.91184756344828699, -0.41052895273466672])
const zcorsika = whitepaper_normal_vec.proj
const xyzcorsika = inv([
  xcorsika.x xcorsika.y xcorsika.z;
  ycorsika.x ycorsika.y ycorsika.z;
  zcorsika.x zcorsika.y zcorsika.z;
 ])

pavel_sim = jldopen("/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.jld2")
config = SimulationConfig(; pavel_sim["config"]...)
injector = Tambo.Injector(config)
geo = Tambo.Geometry(config)
plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)

# This is the parameter we want to vary
ℓ = parse(Float64, ARGS[1]) * units.m

# This is the spacing between detectors
Δs = parse(Float64, ARGS[2]) * units.m # meters

detection_modules = Tambo.make_trianglearray(-2000units.m, 3000units.m, -ℓ/2, ℓ/2, Δs, ϕ=whitepaper_normal_vec.ϕ)
# Remove detectors that are out of range
mask = zmin .< Tambo.plane_z.(getfield.(detection_modules, :x), getfield.(detection_modules, :y), Ref(plane)) .< zmax;
detection_modules = detection_modules[mask]
println(sum(mask))

outdir = ARGS[3]

neutrino_events = np.load("/n/home12/jlazar/CORSIKA_runs.npy")
basedir = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/whitepaper_sim/"
savefile = "$(outdir)/test_events_$(ℓ/units.m)_$(Δs / units.m).npy"

function loadcorsika(files)
    events = Tambo.CorsikaEvent[]
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
      new_poss = []
      for (x, y, z) in zip(xs, ys, zs)
        push!(new_poss, xyzcorsika * [x,y,z])
      end
      xs = getindex.(new_poss, 1)
      ys = getindex.(new_poss, 2)
      zs = getindex.(new_poss, 3)
      ts = x["time"].to_numpy() .* units.second
      ws = x["weight"].to_numpy() .* 1.0
      ids = x["pdg"].to_numpy()
      es = x["kinetic_energy"].to_numpy() .* units.GeV
      for tup in zip(ids, es, xs, ys, zs, ts, ws)
        push!(events, Tambo.CorsikaEvent(tup...))
      end
    end
    return events
end

function did_trigger(files, detection_modules; thresh=1units.m, module_trigger=3, trigger_n=30)
  events = loadcorsika(files)
  # You can't trigger if there are too few particles
  if length(events) <= 1
    return false
  end
  if sum(getfield.(events, :weight)) < trigger_n * 0.5
    return false
  end
  m = zmin .<= getfield.(events, :z) .<= zmax
  hits = Tambo.find_near_hits(
    events,
    detection_modules;
    thresh=thresh
  )
  d = Tambo.make_hit_dict(hits)
  n = 0
  for (_, v) in d
    # a module is not triggered if it has less than three hits
    if v < module_trigger
      continue
    end
    n += v
  end
  return n >= trigger_n
end

thresh = 1units.m # distance from detection module center to consider a hit
triggered_evts = []
for (kdx, neutrino_event) in enumerate(neutrino_events)
  # Find all showers for this events
  files = Tambo.find_extant_files(neutrino_event, basedir)
  if did_trigger(files, detection_modules)
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
