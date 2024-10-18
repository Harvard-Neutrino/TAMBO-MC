module Tambo

export Simulation,
       InjectionConfig,
       #ProposalConfig,
       Geometry,
       CorsikaMap,
       Coord,
       save_simulation, 
       simulator_from_file,
       units,
       coords,
       normal_vecs,
       inject_ν!,
       propagate_τ!,
       #minesite_coord,
       #whitepaper_normal_vec,
       #whitepaper_coord,
       #testsite_coord,
       #minesite_normal_vec,
       #larger_valley_coord,
       #larger_valley_vec,
       inside,
       should_do_corsika,
       latlong_to_xy,
       xy_to_latlong,
       oneweight

using CoordinateTransformations: Translation, AffineMap, LinearMap
using Dierckx: Spline2D
using Distributions: Uniform, Poisson
using JLD2: jldopen, JLDFile
# TODO move h5 to jld2
using HDF5: h5open
using LinearAlgebra: norm
using ProgressBars
using PyCall: PyCall, PyNULL, PyObject
using Random: seed!
using Roots: find_zeros, find_zero
using Rotations: RotX, RotZ
using StaticArrays: SVector, SMatrix 
using TOML

include("units.jl")
include("samplers/angularsamplers.jl")
include("samplers/crosssections.jl")
include("samplers/injectionvolumes.jl")
include("samplers/powerlaws.jl")
include("directions.jl")
include("particles.jl")
include("locations.jl")
include("geometries.jl")
include("tracks.jl")
include("inject.jl")
include("proposal.jl")
include("weightings.jl")
include("taurunner.jl")
include("detector.jl")
include("corsika.jl")

@Base.kwdef mutable struct Simulation
    config::Dict{String, Any}
    results::Dict{String, Any}
end

function relativize!(d::Dict)
    for (k, v) in pairs(d)
        if isa(v, String)
            d[k] = replace(v, "_TAMBO_PATH_" => dirname(pathof(Tambo)))
        elseif isa(v, Dict)
            relativize!(v)
        end
    end
end

function Simulation(config_file::String)
    config = TOML.parsefile(config_file)
    relativize!(config)
    results = Dict{String, Any}()
    return Simulation(config, results)
end

function inject_ν!(sim::Simulation; outkey="injection_events", track_progress=true)
    geo = Geometry(sim.config["geometry"])
    injector = Injector(sim.config["injection"], geo)
    events = Vector{InjectionEvent}(undef, sim.config["steering"]["nevent"])
    itr = 1:sim.config["steering"]["nevent"]
    if track_progress
        itr = ProgressBar(itr)
    end
    for idx in itr
        tr_seed = sim.config["steering"]["seed"] + idx
        event = inject_event(injector, tr_seed)
        events[idx] = event
    end
    sim.results[outkey] = events
end

function propagate_τ!(
    sim::Simulation;
    inkey="injection_events",
    outkey="proposal_events",
    track_progress=true
)
    geo = Geometry(sim.config["geometry"])
    events = Vector{ProposalResult}(undef, sim.config["steering"]["nevent"])
    propagator = ProposalPropagator(sim.config["proposal"])
    injected_events = sim.results[inkey]
    if track_progress
        injected_events = ProgressBar(injected_events)
    end
    for (idx, injected_event) in enumerate(injected_events)
        event = propagator(
            injected_event.final_state,
            geo,
            sim.config["steering"]["seed"] + idx
        )
        events[idx] = event
    end
    sim.results[outkey] = events
end

#function CORSIKAConfig(s::Simulation)
#    propdict = Dict(
#        fn => getfield(s, fn) 
#        for fn in intersect(fieldnames(Simulation), fieldnames(CORSIKAConfig))
#    )
#    return CORSIKAConfig(; propdict...)
#end

function run_airshower!(sim::Simulation; inkey="proposal_events", track_progress=true)
    proposal_events = sim.results[inkey]

    geo = Geometry(sim.config["geometry"])
    # TODO wrap this into a neat little constructor
    plane = Tambo.Plane(
        geo.tambo_normal,
        geo.tambo_coordinates,
        geo
    )
    indices = []
    for (proposal_idx, proposal_event) in enumerate(proposal_events)
        if ~should_do_corsika(proposal_event, plane,geo)
            continue
        end
        for (decay_idx,decay_event) in enumerate(proposal_event.decay_products)
            #wanted to keep indices lined up so checking one at at ime
            if check_neutrino(decay_event)
                continue 
            end 

            push!(indices, [proposal_idx, decay_idx])

            if sim.config["corsika"]["parallelize_corsika"]
                continue 
            end
            corsika_run(
                decay_event,
                sim.config["corsika"],
                geo,
                proposal_idx,
                decay_idx;
                parallelize_corsika=false
            )
        end
    end
    
    if sim.config["corsika"]["parallelize_corsika"]
        corsika_parallel(
            proposal_events,
            geo,
            sim.config["corsika"],
            indices
        )
    end 
    return indices 
end 

function (s::Simulation)(; track_progress=true, should_run_corsika=false)
    throw("Not implemented yet")
end

function save_simulation(s::Simulation, path::String)
    throw("Not implemented yet")
end

end # module
