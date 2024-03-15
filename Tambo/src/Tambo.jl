module Tambo

export SimulationConfig,
       InjectionConfig,
       ProposalConfig,
       Geometry,
       CorsikaMap,
       save_simulation, 
       simulator_from_file,
       units,
       minesite_coord,
       whitepaper_normal_vec,
       whitepaper_coord,
       testsite_coord,
       minesite_normal_vec,
       inside,
       should_do_corsika,
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
include("corsika.jl")
include("weightings.jl")
include("taurunner.jl")
include("detector_responses.jl")

@Base.kwdef mutable struct SimulationConfig
    # General configuration
    n::Int = 10
    seed::Int64 = 925
    run_n::Int64 = 853
    # Geometry configuration
    geo_spline_path::String = realpath("$(@__DIR__)/../../resources/tambo_spline.jld2")
    tambo_coordinates::Coord = minesite_coord
    # Injection configuration
    ν_pdg::Int = 16
    γ::Float64 = 1
    emin::Float64 = 1e6units.GeV
    emax::Float64 = 1e9units.GeV
    θmin::Float64 = 0.0
    θmax::Float64 = π
    ϕmin::Float64 = 0.0
    ϕmax::Float64 = 2π
    r_injection::Float64 = 2000units.m
    l_endcap::Float64 = 1units.km
    diff_xs_path::String = realpath(
        "$(@__DIR__)/../../resources/cross_sections/tables/csms_differential_cdfs.h5"
    )
    # PROPOSAL configuration
    ecut::Float64 = Inf * units.GeV
    vcut::Float64 = 1e-2
    do_interpolate::Bool = true
    do_continuous::Bool = true
    tablespath::String = realpath(
        "$(@__DIR__)/../..//resources/proposal_tables/"
    )

    # CORSIKA configuration 
    parallelize::Bool = true 
    thinning::Float64 = 1e-6 
    hadron_ecut::Float64 = 0.05units.GeV
    em_ecut::Float64 = 0.01units.GeV
    photon_ecut::Float64 = 0.002units.GeV
    mu_ecut::Float64 = 0.05units.GeV 
    shower_dir::String = "showers/"

    injected_events::Vector{InjectionEvent} = InjectionEvent[]
    proposal_events::Vector{ProposalResult} = ProposalResult[]
end

function SimulationConfig(fname::String)
    s = nothing
    jldopen(fname, "r") do f
        s = SimulationConfig(; f["config"]...)
        s.injected_events = f["injected_events"]
        s.proposal_events = f["proposal_events"]
    end
    return s
end

function Injector(config::SimulationConfig)
    geo = Geometry(config)
    cfg = InjectionConfig(config)
    return Injector(cfg, geo)
end

function inject(simulator::SimulationConfig; track_progress=true)
    injector = InjectionConfig(simulator)
    geo = Geometry(simulator)
    return inject(injector, geo, track_progress=track_progress)
end

function ProposalConfig(s::SimulationConfig)
    propdict = Dict(
        fn => getfield(s, fn) 
        for fn in intersect(fieldnames(SimulationConfig), fieldnames(ProposalConfig))
    )
    return ProposalConfig(; propdict...)
end

function CORSIKAConfig(s::SimulationConfig)
    propdict = Dict(
        fn => getfield(s, fn) 
        for fn in intersect(fieldnames(SimulationConfig), fieldnames(CORSIKAConfig))
    )
    return CORSIKAConfig(; propdict...)
end

function InjectionConfig(config::SimulationConfig)
    d = Dict(f=>getfield(config, f) for f in fieldnames(SimulationConfig) if f in fieldnames(InjectionConfig))
    return InjectionConfig(;d...)
end

function Geometry(s::SimulationConfig)
    geo = Geometry(
        s.geo_spline_path,
        s.tambo_coordinates
    )
    #f = jldopen(s.geo_spline_path)
    #spl = f["spline"]
    #mincoord = f["mincoord"]
    #close(f)
    #tambo_xy = latlong_to_xy(s.tambo_coordinates, mincoord)
    #return Geometry(spl, tambo_xy)
    return geo
end

function Base.show(io::IO, s::SimulationConfig)
    print(
        io,
        """
        General configuration
        _____________________
        n: $(s.n)
        seed: $(s.seed)
        run_n: $(s.run_n)

        Geometry configuration
        ______________________
        geo_spline_path: $(s.geo_spline_path)
        tambo_coordinates: $(s.tambo_coordinates)

        Injection configuration
        _______________________
        ν_pdg: $(s.ν_pdg)
        γ: $(s.γ)
        emin: $(s.emin / units.GeV) GeV
        emax: $(s.emax / units.GeV) GeV
        θmin: $(round(s.θmin * 180 / π, sigdigits=3))°
        θmax: $(round(s.θmax * 180 / π, sigdigits=3))°
        ϕmin: $(round(s.ϕmin * 180 / π, sigdigits=3))°
        ϕmax: $(round(s.ϕmax * 180 / π, sigdigits=3))°
        r_injection: $(s.r_injection / units.m) m
        l_endcap: $(s.l_endcap / units.m) m
        diff_xs_path: $(s.diff_xs_path)

        PROPOSAL configuration
        ______________________
        ecut: $(s.ecut / units.GeV) GeV
        vcut: $(s.vcut)
        do_interpolate: $(s.do_interpolate)
        do_continuous: $(s.do_continuous)
        tablespath: $(s.tablespath)

        CORSIKA configuration 
        _____________________
        thinning: $(s.thinning)
        hadron_ecut: $(s.hadron_ecut/ units.GeV) GeV
        muon_ecut: $(s.muon_ecut/ units.GeV) GeV
        em_ecut: $(s.em_ecut/ units.GeV) GeV
        photon_ecut: $(s.photon_ecut/ units.GeV) GeV
        parallelize: $(s.parallelize)
        shower_dir: $(s.shower_dir)
        """
    )
end

function Base.getindex(s::SimulationConfig, fieldstring::String)
    getfield(s, Symbol(fieldstring))
end

function (s::SimulationConfig)(; track_progress=true, should_run_corsika=false)
    seed!(s.seed)
    if track_progress
        println("Making geometry")
    end
    geo = Geometry(
        s.geo_spline_path,
        s.tambo_coordinates
    )
    if track_progress
        println("Injecting events")
    end

    injection_config = InjectionConfig(s)
    injector = Injector(injection_config, geo)
    s.injected_events = injector(track_progress=track_progress)
    if track_progress
        println("Propagating charged leptons")
    end
    proposal_config = ProposalConfig(s)
    propagator = ProposalPropagator(proposal_config)
    s.proposal_events = propagator(
        s.injected_events,
        geo,
        track_progress=track_progress
    )
    if should_run_corsika
        if track_progress
            println("Running CORSIKA showers")
        end
        corsika_config = CORSIKAConfig(s)
        corsika_propagator = CORSIKAPropagator(corsika_config)
        s.corsika_showers = corsika_propagator(
            s.proposal_events, 
            track_progress = track_prograss
        )
    end

end

function dump_to_file(s::SimulationConfig, f::JLDFile)
    resultfields = [:injected_events, :proposal_events]
    f["injected_events"] = s.injected_events
    f["proposal_events"] = s.proposal_events
    f["config"] = Dict(
        Dict(
            fn => getfield(s, fn) for fn in fieldnames(SimulationConfig)
            if fn ∉ resultfields
        )
    )
    return
end

function save_simulation(s::SimulationConfig, path::String)
    @assert(length(s.injected_events)==s.n)
    @assert(length(s.proposal_events)==s.n)
    jldopen(path, "w") do f
        dump_to_file(s, f)
    end
end

end # module
