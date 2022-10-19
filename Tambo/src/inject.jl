include("./samplers/Samplers.jl")
using .Samplers

Base.@kwdef mutable struct Injector
    n::Int = 10
    geo_spline_path::String = realpath("$(@__DIR__)/../../resources/tambo_spline.jld2")
    diff_xs_path::String = realpath(
        "$(@__DIR__)/../../resources/cross_sections/tables/csms_differential_cdfs.h5"
    )
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

function Base.show(io::IO, injector::Injector)
    print(
        io,
        """
        n: $(injector.n)
        seed: $(injector.seed)
        ν_pdg: $(injector.ν_pdg)
        geo_spline_path: $(injector.geo_spline_path)
        diff_xs_path: $(injector.diff_xs_path)
        γ: $(injector.γ)
        emin (GeV): $(injector.emin / units.GeV)
        emax (GeV): $(injector.emax / units.GeV)
        θmin (degrees): $(injector.θmin * 180 / π)
        θmax (degrees): $(injector.θmax * 180 / π)
        ϕmin (degrees): $(injector.ϕmin * 180 / π)
        ϕmax (degrees): $(injector.ϕmax * 180 / π)
        r_injection (m): $(injector.r_injection / units.m)
        l_endcap (m): $(injector.l_endcap / units.m)
        """,
    )
end

function (injector::Injector)(; track_progress=true)
    Random.seed!(injector.seed)
    pl = PowerLaw(injector.γ, injector.emin, injector.emax)
    diff_xs = OutgoingCCEnergy(injector.diff_xs_path, injector.ν_pdg)
    anglesampler = UniformAngularSampler(
        injector.θmin, injector.θmax, injector.ϕmin, injector.ϕmax
    )
    injectionvolume = InjectionVolume(injector.r_injection, injector.l_endcap)
    geo = Geometry(injector.geo_spline_path)
    if track_progress
        iter = ProgressBar(1:(injector.n))
    else
        iter = 1:(injector.n)
    end
    return [
        inject_event(injector.ν_pdg, pl, diff_xs, anglesampler, injectionvolume, geo) for
        _ in iter
    ]
end

"""
    sample_interaction_vertex(
        volume::InjectionVolume,
        p_near::SVector{3},
        d::Direction,
        range,
        geo::Geometry
    )

TBW
"""
function sample_interaction_vertex(
    volume::InjectionVolume, p_near::SVector{3}, d::Direction, range, geo::Geometry
)
    # Make track from point of closest approach to point of entry and exitA
    ti = Track(p_near, d, geo.box)
    to = Track(p_near, reverse(d), geo.box)
    # Compute the intersection of each track with the mountain
    rangesi = computeranges(ti, geo)
    rangeso = computeranges(to, geo)
    # Compute the colum depth for both incoming and outgoing portions
    cdi = endcapcolumndepth(ti, volume.l_endcap, range, rangesi)
    cdo = endcapcolumndepth(to, volume.l_endcap, 0.0, rangeso)
    # sample column depth uniformly and subtract incoming column depth
    cd = rand(Uniform(-cdo, cdi))
    # If the remainder is positive, you need to be in incoming track, else outgoing
    cd > 0 ? tr = ti : tr = to
    cd > 0 ? ranges = rangesi : ranges = rangeso
    # Find affine parameter where we have traversed proper column depth
    λ_int = inversecolumndepth(tr, abs(cd), geo, ranges)
    # Convert affine parameter to a physical location
    p_int = tr(λ_int)
    return p_int
end

"""
    endcapcolumndepth(t::Track, l_endcap::Float64, range::Float64, ranges::Vector)

TBW
"""
function endcapcolumndepth(t::Track, l_endcap::Float64, range::Float64, ranges::Vector)
    cd = totalcolumndepth(t, ranges)
    if t.norm <= l_endcap
        cd_endcap = cd
    else
        cd_endcap = minimum([columndepth(t, l_endcap / t.norm, ranges) + range, cd])
    end
    return cd_endcap
end

struct InjectionEvent
    initial_state::Particle
    final_state::Particle
end

function Base.show(io::IO, event::InjectionEvent)
    print(
        io,
        """
        initial_state:
        $(event.initial_state)

        final_state:
        $(event.final_state)
        """,
    )
end

"""
    inject_event(injector::Injector)

TBW
"""
function inject_event(
    ν_pdg::Int,
    e_sampler,
    diff_xs::OutgoingCCEnergy,
    anglesampler,
    injectionvolume::InjectionVolume,
    geo::Geometry,
)
    e_init = rand(e_sampler)
    # Sometimes this sampler runs into floating poiting precision issues that give
    # an energy higher than the initial energy
    # I should probably figure that out but for now I leave it.
    # Hopefully the splines fix this
    e_final = rand(diff_xs, e_init)
    #e_final = minimum([rand(diff_xs, e_init), e_init])
    range = lepton_range(e_init, abs(ν_pdg) == 16)
    d = Direction(rand(anglesampler)...)
    # Construct roatation to plane perpindicular to direction
    r = (Rotations.RotX(d.θ) * RotZ(π / 2 - d.ϕ))'
    p_near = r * rand(injectionvolume)
    p_int = sample_interaction_vertex(injectionvolume, p_near, d, range, geo)
    initial_state = Particle(ν_pdg, e_init, p_int, d, nothing)
    final_state = Particle(ν_pdg - sign(ν_pdg), e_final, p_int, d, initial_state)
    event = InjectionEvent(initial_state, final_state)
    return event
end

### Convenience functions ###

function Base.getindex(v::Vector{InjectionEvent}, s::String)
    return getfield.(v, Symbol(s))
end
