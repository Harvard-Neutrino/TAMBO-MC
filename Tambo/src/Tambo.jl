module Tambo

export Simulation,
       Geometry,
       CorsikaMap,
       Coord,
       units,
       coords,
       normal_vecs,
       inject_ν!,
       propagate_τ!,
       run_airshower!,
       oneweight,
       save_simulation

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
    function Simulation(config, results)
        @assert "geometry" in keys(config) "Geometry information must be provided"
        @assert "steering" in keys(config) "Steering information must be provided"
        return new(config, results)
    end
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

function inject_ν!(
    sim::Simulation,
    config::Dict{String, Any};
    outkey="injection",
    track_progress=true
)
    relativize!(config)
    sim.config[outkey] = config
    geo = Geometry(sim.config["geometry"])
    injector = Injector(config, geo)
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

function inject_ν!(
    sim::Simulation,
    config_file::String;
    outkey="injection",
    track_progress=true
)
    config = relativize!(TOML.parsefile(config_file))
    inject_ν!(sim, config; outkey=outkey, track_progress=track_progress)
end

function propagate_τ!(
    sim::Simulation,
    config::Dict{String, Any};
    inkey="injection",
    outkey="proposal",
    track_progress=true
)
    relativize!(config)
    sim.config[outkey] = config
    geo = Geometry(sim.config["geometry"])
    events = Vector{ProposalResult}(undef, sim.config["steering"]["nevent"])
    propagator = ProposalPropagator(config)
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

function propagate_τ!(
    sim::Simulation,
    config_file::String;
    inkey::String="injection",
    outkey::String="proposal",
    track_progress::Bool=true
)
    config = relativize!(TOML.parsefile(config_file))
    propagate_τ!(sim, config; inkey=inkey, outkey=outkey, track_progress=track_progress)
end

function run_airshower!(
    sim::Simulation,
    config::Dict{String, Any};
    outkey="corsika",
    inkey="proposal",
    track_progress=true
)
    relativize!(config)
    proposal_events = sim.results[inkey]
    
    sim.config[outkey] = config
    geo = Geometry(sim.config["geometry"])
    # TODO wrap this into a neat little constructor
    plane = Tambo.Plane(
        geo.tambo_normal,
        geo.tambo_coordinates,
        geo
    )
    indices = Vector{Tuple{Int64, Int64}}()
    for (proposal_idx, proposal_event) in enumerate(proposal_events)
        if ~should_do_corsika(proposal_event, plane,geo)
            continue
        end
        for (decay_idx,decay_event) in enumerate(proposal_event.decay_products)
            #wanted to keep indices lined up so checking one at at ime
            if check_neutrino(decay_event)
                continue 
            end 

            push!(indices, (proposal_idx, decay_idx))

            if sim.config[outkey]["parallelize_corsika"]
                continue 
            end
            corsika_run(
                decay_event,
                sim.config[outkey],
                geo,
                proposal_idx,
                decay_idx;
                parallelize_corsika=false
            )
        end
    end
    
    if sim.config["corsika"]["parallelize_corsika"]
        println(indices)
        corsika_parallel(
            proposal_events,
            geo,
            sim.config["corsika"],
            indices
        )
    end 
    sim.results[outkey] = indices
    #return indices 
end 

function run_airshower!(sim::Simulation, config_file::String; outkey="corsika", inkey="proposal_events", track_progress=true)
    config = relativize!(TOML.parsefile(config_file))
    run_airshower!(sim, config; outkey=outkey, inkey=inkey, track_progress=track_progress)
end

function (s::Simulation)(; track_progress=true, should_run_corsika=false)
    throw("Not implemented yet")
end

function dump_to_file(s::Simulation, f::JLDFile)
    #resultfields = [:injected_events, :proposal_events, :corsika_indices]
    println("Keys $(keys(s.results))")
    println(s.results)
    f["injected_events"] = s.results["injection"]
    f["proposal_events"] = s.results["proposal"]
    f["corsika_indices"] = s.results["corsika"]
    #f["config"] = Dict(
    #    Dict(
    #        fn => getfield(s, fn) for fn in fieldnames(SimulationConfig)
    #        if fn ∉ resultfields
    #    )
    #)
    f["config"] = s.config
    return
end

function save_simulation(s::Simulation, path::String)
    @assert length(s.results["injection"]) == s.config["steering"]["nevent"]
    @assert length(s.results["proposal"]) == s.config["steering"]["nevent"]
    jldopen(path, "w") do file
        dump_to_file(s, file)
    end
end

end # module