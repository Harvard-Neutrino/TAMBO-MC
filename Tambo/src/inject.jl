Base.@kwdef mutable struct InjectionConfig
    n::Int = 10
    xs_dir::String = realpath(
        "$(@__DIR__)/../../resources/cross_sections/tables/"
    )
    xs_model::String = "csms"
    interaction::Interaction = Interaction(1) # Charged current
    tambo_coordinates::Coord = minesite_coord
    ν_pdg::Int = 16
    γ::Float64 = 1
    emin::Float64 = 1e6units.GeV
    emax::Float64 = 1e9units.GeV
    θmin::Float64 = 0.0
    θmax::Float64 = π
    ϕmin::Float64 = 0.0
    ϕmax::Float64 = 2π
    r_injection::Float64 = 900units.m
    l_endcap::Float64 = 1units.km
    seed::Int64 = 0
end

struct Injector
    config::InjectionConfig
    powerlaw::PowerLaw
    xs::CrossSection
    anglesampler::UniformAngularSampler
    injectionvolume::SymmetricInjectionCylinder
    geo::Geometry
end

function f()
    return 4
end

function Injector(config::InjectionConfig, geo::Geometry)
    pl = PowerLaw(config.γ, config.emin, config.emax)
    xs = CrossSection(config.xs_dir, config.xs_model, config.ν_pdg, config.interaction)
    anglesampler = UniformAngularSampler(
        config.θmin,
        config.θmax,
        config.ϕmin,
        config.ϕmax
    )

    injectionvolume = SymmetricInjectionCylinder(config.r_injection, config.l_endcap)
    return Injector(config, pl, xs, anglesampler, injectionvolume, geo)
end

function inject_event(injector::Injector)
    event = inject_event(
        injector.config.ν_pdg,
        injector.powerlaw,
        injector.xs,
        injector.anglesampler,
        injector.injectionvolume,
        injector.geo
    )
    return event
end

function (injector::Injector)(; track_progress=true)
    seed!(injector.config.seed)
    iter = 1:injector.config.n
    if track_progress
        iter = ProgressBar(iter)
    end
    return [inject_event(injector) for _ in iter]
end

struct InjectionEvent
    entry_state::Particle
    initial_state::Particle
    final_state::Particle
    physX::Float64
    genX::Float64
end

"""
    sample_interaction_vertex(
        volume::InjectionVolume,
        closest_approach::SVector{3},
        d::Direction,
        range,
        geo::Geometry
    )

TBW
"""
function sample_interaction_vertex(
    volume::SymmetricInjectionCylinder,
    closest_approach::SVector{3},
    d::Direction,
    range::Float64,
    geo::Geometry
)
    track = Track(closest_approach, reverse(d), geo.box)
    segments = computesegments(track, geo)
    tot_X = endcapcolumndepth(track, volume.l_endcap, range, segments)
    X = rand(Uniform(0.0, tot_X))
    λ_int = inversecolumndepth(track, X, geo, segments)
    p_int = track(λ_int)
    return p_int, tot_X
end

"""
    endcapcolumndepth(t::Track, l_endcap::Float64, range::Float64, segments::Vector{Segment})

TBW
"""
function endcapcolumndepth(t::Track, l_endcap::Float64, range::Float64, segments::Vector{Segment})
    cd = totalcolumndepth(t, segments)
    if t.norm <= l_endcap
        return cd
    end
    cd_endcap = minimum([columndepth(t, l_endcap / t.norm, segments) + range, cd])
    return cd_endcap
end

"""
    inject_event(
    ν_pdg::Int,
    e_sampler,
    diff_xs::CrossSection,
    anglesampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry,
)
    ν_pdg::Int,
    e_sampler,
    diff_xs::OutgoingEnergy,
    anglesampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry,
)

TBW
"""
function inject_event(
    ν_pdg::Int,
    power_law::PowerLaw,
    xs::CrossSection,
    anglesampler::UniformAngularSampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry,
)
    direction = Direction(rand(anglesampler)...)
    # Rotation to plane perpindicular to direction
    rotator = (RotX(direction.θ) * RotZ(π / 2 - direction.ϕ))'
    closest_approach = rotator * rand(injectionvolume)
    xb = intersect(closest_approach, reverse(direction), geo.box)

    proposed_e_init = rand(power_law)
    proposed_particle = Particle(ν_pdg, proposed_e_init, xb, direction, nothing)
    particle_entry, physX = tr_propagate(proposed_particle, geo.tambo_offset.z)
    e_final = rand(xs, particle_entry.energy)

    range = lepton_range(e_final, ν_pdg - sign(ν_pdg))
    p_int, genX = sample_interaction_vertex(injectionvolume, closest_approach, direction, range, geo)
    final_state = Particle(ν_pdg - sign(ν_pdg), e_final, p_int, direction, particle_entry)
    # Now our generation column depth is the same as our physical column depth so we set them equal
    event = InjectionEvent(
        particle_entry,
        proposed_particle,
        final_state,
        genX,
        genX
    )
    return event
end

function save_simulation(config::InjectionConfig, path::String)
    @assert length(s.injected_events)==s.n "Looks like you didn't inject the right number of events"
    jldopen(path, "w") do f
        dump_to_file(s, f)
    end
end

### Convenience functions ###

function Base.getindex(v::Vector{InjectionEvent}, s::String)
    return getfield.(v, Symbol(s))
end

function Base.show(io::IO, config::InjectionConfig)
    print(
        io,
        """
        n=$(config.n)
        seed=$(config.seed)
        ν_pdg=$(config.ν_pdg)
        interaction=$(config.interaction)
        xs_model=$(config.xs_model)
        γ=$(config.γ)
        emin=$(config.emin / units.GeV) GeV
        emax=$(config.emax / units.GeV) GeV
        θmin=$(round(config.θmin * 180 / π, sigdigits=3))°
        θmax=$(round(config.θmax * 180 / π, sigdigits=3))°
        ϕmin=$(round(config.ϕmin * 180 / π, sigdigits=3))°
        ϕmax=$(round(config.ϕmax * 180 / π, sigdigits=3))°
        r_injection=$(config.r_injection / units.m) m 
        l_endcap=$(config.l_endcap / units.m) m
        xs_dir=$(config.xs_dir)
        """
    )
end
