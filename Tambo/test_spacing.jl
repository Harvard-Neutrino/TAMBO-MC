using Revise
using Pkg
Pkg.activate(".")
using Tambo
using PyCall
using JLD2
using CSV

np = pyimport("numpy")
ak = pyimport("awkward")

pavel_sim = jldopen("/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.jld2")
config = SimulationConfig(; pavel_sim["config"]...)
injector = Tambo.Injector(config)
geo = Tambo.Geometry(config)
plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)

# This rotates the geometry to align with the river
#ϕ = atan(-0.366163, whitepaper_normal_vec.)
ϕ = atan(-0.366163, 0.452174)

# This is the parameter we want to vary
ℓ = parse(Float64, ARGS[1]) * units.m

# This is the spacing between detectors
Δs = parse(Float64, ARGS[2]) * units.m # meters

detection_modules = Tambo.make_trianglearray(-800units.m, 2000units.m, -ℓ/2, ℓ/2, Δs, ϕ=ϕ)

# Remove detectors that are out of range
zmin = 2100units.m - geo.tambo_offset.z
zmax = 5000units.m - geo.tambo_offset.z
mask = zmin .< Tambo.plane_z.(getfield.(detection_modules, :x), getfield.(detection_modules, :y), Ref(plane), Ref(geo)) .< zmax
detection_modules = detection_modules[mask]

outdir = ARGS[3]

neutrino_events = np.load("/n/home12/jlazar/CORSIKA_runs.npy")
basedir = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/whitepaper_sim/"
pavel_magic_csv = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.csv"
savefile = "$(outdir)/test_events_$(ℓ/units.m)_$(Δs / units.m).npy"

csv = CSV.File(pavel_magic_csv)
idxs = getindex.(csv, 16)
subidxs = getindex.(csv, 17)

thresh = 1units.m # distance from detection module center to consider a hit
triggered_evts = []
for (kdx, neutrino_event) in enumerate(neutrino_events)
  # Find all showers for this events
  files = Tambo.find_extant_files(neutrino_event, basedir)
  # Make a vactor that we will fill with hits
  plz = Tambo.Hit[]
  for file in files
    desc = split(file, "/")[end-2]
    idx = parse(Int, split(desc, "_")[end-1])
    subidx = parse(Int, split(desc, "_")[end])
    row = first(csv[(idx.==idxs) .&& (subidx.==subidxs)])

    # try to load the file
    # We need this try/catch logic cuz empty files give issues
    a = nothing
    try
      a = ak.from_parquet(file)
    catch
      continue
    end
    # Find hits that fall within a certain distance of any detection unit
    hits = Tambo.find_near_hits(
      a,
      detection_modules,
      row.inter_x*units.km,
      row.inter_y * units.km;
      thresh=thresh
    )
    # Add these hits to the global list of such hits
    # plz = hcat(plz, hits) Is probably right
    for hit in hits
      push!(plz, hit)
    end
  end
  # A dictionary whose keys are detection modules and whose values are
  # integer number of particles. This does Poisson sampling for events
  # that have thinning
  d = Tambo.make_hit_dict(plz)
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
  if n >= 30
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
