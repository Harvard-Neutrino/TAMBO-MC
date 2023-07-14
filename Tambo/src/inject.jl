using .Samplers

Base.@kwdef mutable struct InjectionConfig
    n::Int = 10
    diff_xs_path::String = realpath(
        "$(@__DIR__)/../../resources/cross_sections/tables/csms_differential_cdfs.h5"
    )
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

function Base.show(io::IO, injector::InjectionConfig)
    print(
        io,
        """
        n: $(injector.n)
        seed: $(injector.seed)
        ν_pdg: $(injector.ν_pdg)
        diff_xs_path: $(injector.diff_xs_path)
        γ: $(injector.γ)
        emin (GeV): $(injector.emin / units.GeV)
        emax (GeV): $(injector.emax / units.GeV)
        θmin (degrees): $(round(injector.θmin * 180 / π, sigdigits=3))°
        θmax (degrees): $(round(injector.θmax * 180 / π, sigdigits=3))°
        ϕmin (degrees): $(round(injector.ϕmin * 180 / π, sigdigits=3))°
        ϕmax (degrees): $(round(injector.ϕmax * 180 / π, sigdigits=3))°
        r_injection (m): $(injector.r_injection / units.m)
        l_endcap (m): $(injector.l_endcap / units.m)"""
    )
end

function inject(injector::InjectionConfig, geo::Geometry; track_progress=true)
    seed!(injector.seed)
    spectrum = PowerLaw(injector.γ, injector.emin, injector.emax)
    diff_xs = OutgoingCCEnergy(injector.diff_xs_path, injector.ν_pdg)
    anglesampler = UniformAngularSampler(
        injector.θmin,
        injector.θmax,
        injector.ϕmin,
        injector.ϕmax
    )

    injectionvolume = SymmetricInjectionCylinder(injector.r_injection, injector.l_endcap)

    iter = 1:(injector.n)
    if track_progress
        iter = ProgressBar(iter)
    end

    return [
        inject_event(injector.ν_pdg, spectrum, diff_xs, anglesampler, injectionvolume, geo) for
        _ in iter
    ]
end

function (injector::InjectionConfig)(geo; track_progress=true)
    return inject(injector, geo, track_progress=track_progress)
end

struct InjectionEvent
    entry_state::Particle
    initial_state::Particle
    final_state::Particle
    X::Float64
end

function Base.show(io::IO, event::InjectionEvent)
    print(
        io,
        """
        initial_state:
        $(event.initial_state)

        final_state:
        $(event.final_state)""",
    )
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
    # Track from closest approach to incoming edge
    track = Track(closest_approach, d, geo.box)
    # Only computing these once speeds things up
    segments = computesegments(track, geo)
    tot_X = endcapcolumndepth(track, volume.l_endcap, range, segments)
    X = rand(Uniform(0.0, tot_X))
    λ_int = inversecolumndepth(track, X, geo, segments)
    p_int = track(λ_int)
    return p_int
end

"""
    endcapcolumndepth(t::Track, l_endcap::Float64, range::Float64, segments::Vector{Segment})

TBW
"""
function endcapcolumndepth(t::Track, l_endcap::Float64, range::Float64, segments::Vector{Segment})
    cd = totalcolumndepth(t, segments)
    if t.norm <= l_endcap
        cd_endcap = cd
    else
        cd_endcap = minimum([columndepth(t, l_endcap / t.norm, segments) + range, cd])
    end
    return cd_endcap
end

"""
    inject_event(
    ν_pdg::Int,
    e_sampler,
    diff_xs::OutgoingCCEnergy,
    anglesampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry,
)

TBW
"""
function inject_event(
    ν_pdg::Int,
    e_sampler,
    diff_xs::OutgoingCCEnergy,
    anglesampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry,
)
    direction = Direction(rand(anglesampler)...)
    # Construct roatation to plane perpindicular to direction
    rotator = (RotX(direction.θ) * RotZ(π / 2 - direction.ϕ))'
    closest_approach = rotator * rand(injectionvolume)
    xb = intersect(closest_approach, reverse(direction), geo.box)

    proposed_e_init = rand(e_sampler)
    proposed_particle = Particle(ν_pdg, proposed_e_init, xb, direction, nothing)
    particle_entry, X = tr_propagate(proposed_particle, geo.tambo_offset.z)
    e_final = rand(diff_xs, particle_entry.energy)

    range = lepton_range(particle_entry.energy, ν_pdg)
    p_int = sample_interaction_vertex(injectionvolume, closest_approach, direction, range, geo)
    final_state = Particle(ν_pdg - sign(ν_pdg), e_final, p_int, direction, particle_entry)
    event = InjectionEvent(proposed_particle, particle_entry, final_state, X)
    return event
end

function save_simulation(injector::InjectionConfig, path::String)
    @assert length(s.injected_events)==s.n "Looks like you didn't inject the right number of events"
    jldopen(path, "w") do f
        dump_to_file(s, f)
    end
end

### Convenience functions ###

function Base.getindex(v::Vector{InjectionEvent}, s::String)
    return getfield.(v, Symbol(s))
end
