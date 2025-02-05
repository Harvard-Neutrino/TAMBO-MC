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
       identify_taus_to_shower!,
       shower_taus!,
       run_subshower!,
       run_airshower!,
       oneweight,
       save_simulation

using CoordinateTransformations: Translation, AffineMap, LinearMap
using Dierckx: Spline2D
using Distributions: Uniform, Poisson
using JLD2: jldopen, JLDFile, load
# TODO move h5 to jld2
using HDF5: h5open
using LinearAlgebra: norm
using ProgressBars
using PyCall: PyCall, PyNULL, PyObject
using Random: seed!, rand
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
            d[k] = replace(v, "_TAMBOSIM_PATH_" => ENV["TAMBOSIM_PATH"])
        elseif isa(v, Dict)
            relativize!(v)
        end
    end
end

function validate_config_file(config::Dict{String, Any})
    # Check that only expected configuration parameters are present
    # so user doesn't think they're setting parameters they aren't

    expected_top_level_keys = Set(["geometry", "steering", "injection", "proposal", "corsika"])
    unexpected_keys = setdiff(Set(keys(config)), expected_top_level_keys)
    if !isempty(unexpected_keys)
        error("Unexpected keys found in config file: ", unexpected_keys)
    end

    expected_steering_keys = Set(["nevent", "pinecone", "run_number"])
    unexpected_keys = setdiff(Set(keys(config["steering"])), expected_steering_keys)
    if !isempty(unexpected_keys)
        error("Unexpected keys found in steering section of config file: ", unexpected_keys)
    end

    expected_geo_keys = Set(["geo_spline_path", "tambo_coordinates", "plane_orientation"])
    unexpected_keys = setdiff(Set(keys(config["geometry"])), expected_geo_keys)
    if !isempty(unexpected_keys)
        error("Unexpected keys found in geometry section of config file: ", unexpected_keys)
    end

    expected_injection_keys = Set(["nu_pdg", "gamma", "gamma", "emin", "emax", "thetamin", "thetamax", "phimin", "phimax", "r_injection", "l_endcap", "xs_dir", "xs_model", "interaction", "track_progress"])
    unexpected_keys = setdiff(Set(keys(config["injection"])), expected_injection_keys)
    if !isempty(unexpected_keys)
        error("Unexpected keys found in injection section of config file: ", unexpected_keys)
    end

    expected_proposal_keys = Set(["ecut", "vcut", "do_interpolate", "do_continuous", "tablespath", "track_progress"])
    unexpected_keys = setdiff(Set(keys(config["proposal"])), expected_proposal_keys)
    if !isempty(unexpected_keys)
        error("Unexpected keys found in proposal section of config file: ", unexpected_keys)
    end

    expected_corsika_keys = Set(["should_run_corsika", "parallelize_corsika", "thinning", "hadron_ecut", "em_ecut", "photon_ecut", "mu_ecut", "corsika_path", "track_progress", "FLUPRO", "FLUFOR"])
    unexpected_keys = setdiff(Set(keys(config["corsika"])), expected_corsika_keys)
    if !isempty(unexpected_keys)
        error("Unexpected keys found in corsika section of config file: ", unexpected_keys)
    end
end

function Simulation(config_file::String, injection_file::String="")
    config = TOML.parsefile(config_file)
    validate_config_file(config)
    relativize!(config)
    if injection_file == ""
        results = Dict{String, Any}()

    else
        results = load(injection_file)
    end
    return Simulation(config, results)
end

function inject_ν!(
    sim::Simulation,
    config::Dict{String, Any},
    simset_id::String,
    seed::Int64;
    outkey="injected_events",
    track_progress=true
)
    relativize!(config)
    sim.config[outkey] = config

    seed!(seed)
    track_progress = sim.config[outkey]["track_progress"]

    geo = Geometry(sim.config["geometry"])
    injector = Injector(config, geo)
    events = Vector{InjectionEvent}(undef, sim.config["steering"]["nevent"])
    itr = 1:sim.config["steering"]["nevent"]
    if track_progress
        println("Injecting neutrinos")
        itr = ProgressBar(itr)
    end

    event_id_offset = parse(simset_id, Int) * sim.config["steering"]["nevent"]
    for idx in itr
        tr_seed = seed + idx
        event_id = event_id_offset + idx
        event = inject_event(injector, event_id, tr_seed)
        events[idx] = event
    end
    sim.results[outkey] = events
end

function inject_ν!(
    sim::Simulation,
    config_file::String;
    outkey="injected_events",
    track_progress=true
)
    config = relativize!(TOML.parsefile(config_file))
    inject_ν!(sim, config; outkey=outkey, track_progress=track_progress)
end

function propagate_τ!(
    sim::Simulation,
    config::Dict{String, Any},
    seed::Int64;
    inkey="injected_events",
    outkey="proposal_events",
    track_progress=true
)
    relativize!(config)
    sim.config[outkey] = config

    seed!(seed)
    track_progress = sim.config[outkey]["track_progress"]

    geo = Geometry(sim.config["geometry"])
    events = Vector{ProposalResult}(undef, sim.config["steering"]["nevent"])
    propagator = ProposalPropagator(config)
    injected_events = sim.results[inkey]
    if track_progress
        println("Propagating taus")
        injected_events = ProgressBar(injected_events)
    end
    for (idx, injected_event) in enumerate(injected_events)
        event = propagator(
            injected_event.final_state,
            geo,
            round(Int, rand()) + idx
        )
        events[idx] = event
    end
    sim.results[outkey] = events
end

function propagate_τ!(
    sim::Simulation,
    config_file::String;
    inkey::String="injected_events",
    outkey::String="proposal_events",
    track_progress::Bool=true
)
    config = relativize!(TOML.parsefile(config_file))
    propagate_τ!(sim, config; inkey=inkey, outkey=outkey, track_progress=track_progress)
end

function identify_taus_to_shower!(
    sim::Simulation,
    config::Dict{String, Any},
    seed::Int64;
    inkey="proposal_events",
    outkey="corsika_indices",
    track_progress=true
)
    relativize!(config)
    proposal_events = sim.results[inkey]

    track_progress = sim.config["corsika"]["track_progress"]
    
    sim.config[outkey] = config
    geo = Geometry(sim.config["geometry"])
    # TODO wrap this into a neat little constructor
    plane = Tambo.Plane(
        geo.tambo_normal,
        geo.tambo_coordinates,
        geo
    )

    indices = Vector{Tuple{Int64, Int64}}()

    if track_progress
        println("Identifying taus to shower")
        proposal_events = ProgressBar(proposal_events)
    end
    for (proposal_idx, proposal_event) in enumerate(proposal_events)
        if ~should_do_corsika(proposal_event, plane,geo)
            continue
        end
        for (decay_idx,decay_event) in enumerate(proposal_event.decay_products)
            #wanted to keep indices lined up so checking one at at ime
            if check_neutrino(decay_event)
                continue 
            end
            event_id = proposal_idx * sim.config["steering"]["nevent"] 
            push!(indices, (event_id, decay_idx))
        end
    end
    sim.results[outkey] = indices
end

function shower_taus!(
    sim::Simulation,
    config::Dict{String, Any};
    proposal_ids_key="corsika_indices",
    proposal_events_key="proposal_events",
    track_progress=true
)
    relativize!(config)

    # TODO: I think we should seed in the script calling this function,
    # not in the function itself
    seed!(sim.config["steering"]["seed"])

    # Get proposal event ids corrosponding to the taus that passed should_do_corsika
    should_do_corsika_proposal_ids = sort(unique([t[1] for t in sim.results[proposal_ids_key]]))
    proposal_events = sim.results[proposal_events_key]
    
    sim.config["corsika"] = config
    geo = Geometry(sim.config["geometry"])
    # TODO wrap this into a neat little constructor
    plane = Tambo.Plane(
        geo.tambo_normal,
        geo.tambo_coordinates,
        geo
    )
    indices = Vector{Tuple{Int64, Int64}}()
    # TODO: this double for loop is duplicating work done in identify_taus_to_shower!
    for (proposal_idx, proposal_event) in zip(should_do_corsika_proposal_ids, proposal_events[should_do_corsika_proposal_ids])
        for (decay_idx,decay_event) in enumerate(proposal_event.decay_products)
            #wanted to keep indices lined up so checking one at at ime
            push!(indices, (proposal_idx, decay_idx))

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
        println(indices)
        corsika_parallel(
            proposal_events,
            geo,
            sim.config["corsika"],
            indices
        )
    end 
end 

function run_subshower!(
    sim::Simulation,
    config::Dict{String, Any},
    proposal_id::Int64,
    decay_id::Int64,
    seed::Int64;
    #output_path::String;
    proposal_ids_key="corsika_indices",
    proposal_events_key="proposal_events",
    track_progress=true
)
    relativize!(config)

    seed!(seed)

    sim.config["corsika"] = config
    geo = Geometry(sim.config["geometry"])
    plane = Tambo.Plane(
        geo.tambo_normal,
        geo.tambo_coordinates,
        geo
    )

    corsika_run(
        sim.results[proposal_events_key][proposal_id].decay_products[decay_id],
        sim.config["corsika"],
        geo,
        proposal_id,
        decay_id,
        seed,
        parallelize_corsika=false
    )
end

function run_airshower!( # TODO: obsolete?
    sim::Simulation,
    config::Dict{String, Any};
    outkey="corsika_indices",
    inkey="proposal_events",
    track_progress=true
)
    relativize!(config)
    proposal_events = sim.results[inkey]

    # TODO: I think we should seed in the script calling this function,
    # not in the function itself
    seed!(sim.config["steering"]["seed"])
    
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
end 

function run_airshower!(sim::Simulation, config_file::String; outkey="corsika_indices", inkey="proposal_events", track_progress=true)
    config = relativize!(TOML.parsefile(config_file))
    run_airshower!(sim, config; outkey=outkey, inkey=inkey, track_progress=track_progress)
end

function (s::Simulation)(; track_progress=true, should_run_corsika=false)
    # TODO: I think we should seed in the script calling this function,
    # not in the function itself
    seed!(s.config["steering"]["seed"])

    if track_progress
        println("Constructing geometry")
    end
    geo = Geometry(s.config["geometry"])

    if track_progress
        println("Injecting neutrinos")
    end
    inject_ν!(s, s.config["injected"], track_progress=track_progress)

    if track_progress
        println("Propagating taus")
    end
    propagate_τ!(s, s.config["proposal"], track_progress=track_progress)

    if should_run_corsika
        if track_progress
            println("Running airshowers")
        end
        run_airshower!(s, s.config["corsika"], track_progress=track_progress)
    end
end

function dump_to_file(s::Simulation, f::JLDFile)
    #resultfields = [:injected_events, :proposal_events, :corsika_indices]
    f["injected_events"] = s.results["injected_events"]
    f["proposal_events"] = s.results["proposal_events"]
    f["corsika_indices"] = s.results["corsika_indices"]
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
    @assert length(s.results["injected_events"]) == s.config["steering"]["nevent"]
    @assert length(s.results["proposal_events"]) == s.config["steering"]["nevent"]
    jldopen(path, "w") do file
        dump_to_file(s, file)
    end
end

end # module