mutable struct Injector
    n::Int
    geo_spline_path::String
    diff_xs_path::String
    ν_pdg::Int
    γ::Float64
    emin::Float64
    emax::Float64
    θmin::Float64
    θmax::Float64
    ϕmin::Float64
    ϕmax::Float64
    r_injection::Float64
    l_endcap::Float64
    seed::Int64
end

function Injector()
    n = 10
    geo_spline_path = realpath("$(@__DIR__)/../../resources/tambo_spline.jld2")
    xs_path = realpath("$(@__DIR__)/../../resources/cross_sections/tables/csms_differential.h5")
    ν_pdg = 16
    γ = 1
    emin = 1e6units[:GeV]
    emax = 1e9units[:GeV]
    θmin = 0
    θmax = π
    ϕmin = 0
    ϕmax = 2π
    r_injection = 900units[:m]
    l_endcap = 1units[:km]
    seed = 0
    Injector(
        n,
        geo_spline_path,
        xs_path,
        ν_pdg,
        γ,
        emin,
        emax,
        θmin,
        θmax,
        ϕmin,
        ϕmax,
        r_injection,
        l_endcap,
        seed
    )
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
    θmin: $(injector.θmin)
    θmax: $(injector.θmax)
    ϕmin: $(injector.ϕmin)
    ϕmax: $(injector.ϕmax)
    r_injection (m): $(injector.r_injection / units.m)
    l_endcap (m): $(injector.l_endcap / units.m)
    """
    )
end

function (injector::Injector)(track_progress=true)
    Random.seed!(injector.seed)
    pl = PowerLaw(injector.γ, injector.emin, injector.emax)
    geo = Geometry(injector.geo_spline_path)
    diff_xs = DifferentialXS(injector.diff_xs_path, injector.ν_pdg)
    if track_progress
        iter = ProgressBar(1:injector.n)
    else
        iter = 1:injector.n
    end
    [inject_event(injector, pl, geo, diff_xs) for _ in iter]
end

"""
    perpendicular_plane(θ, ϕ, b, ψ)

rotates the vector in the xy-plane defined by (`b`, `ψ`) to a plane
perpendicular to the 3D unit vector defined by (`θ`, ϕ). `return_transform`
returns the rotation matrix as well as the transformed vector

# Example
```julia-repl
julia> pv = perpendicular_plane(π/3, 5π/4, 200, 7π/6)
3-element SVector{3, Float64} with indices SOneTo(3):
 157.82982619848627
 -87.11914807983158
  86.6025403784438

julia> sum(pv .* [sin(π/3)cos(5π/4), sin(π/3)sin(5π/4), cos(π/3)])
-7.105427357601002e-15
```
"""
function perpendicular_plane(θ, ϕ, b, ψ)
    # Construct vector in the plane of normal coordinate system
    bv = SVector{3}([b*cos(ψ), b*sin(ψ), b*0])
    r = (Rotations.RotX(θ) * RotZ(π/2-ϕ))'
    r * bv
end

function sample_interaction_vertex(
    injector::Injector,
    p_near::SVector{3},
    d::Direction,
    range,
    geo::Geometry,
)
    # Make track from point of closest approach to point of entry and exitA
    ti = Track(p_near, d, geo.box)
    to = Track(p_near, reverse(d), geo.box)
    # Compute the intersection of each track with the mountain
    rangesi = computeranges(ti, geo)
    rangeso = computeranges(to, geo)
    # Compute the colum depth for both incoming and outgoing portions
    # TODO figure out why the same tracks give different cds...
    # Specifically this happens when you choose θmin = θmax = π
    cdi = endcapcolumndepth(ti, injector.l_endcap, range, rangesi)
    cdo = endcapcolumndepth(to, injector.l_endcap, 0.0, rangeso)
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
        cd_endcap = minimum(
            [columndepth(t, l_endcap/t.norm, ranges) + range, cd]
        )
    end
    cd_endcap
end

struct InjectionEvent
    initial_state::Particle
    final_state::Particle
    proposal_result::ProposalResult
end

function Base.show(io::IO, event::InjectionEvent)
    print(
        io,
        """
        initial_state: $(event.initial_state)
        final_state: $(event.final_state)
        proposal_result $(event.proposal_result)
end
        """
    )
end

"""
    inject_event(injector::Injector)

TBW
"""
function inject_event(
    injector::Injector,
    e_sampler,
    geo::Geometry,
    diff_xs::DifferentialXS
)
    # Sample initial and final state energies
    e_init = rand(e_sampler)
    # Sometimes this sampler runs into floating poiting precision issues that give
    # an energy higher than the initial energy
    e_final = maximum([sample(diff_xs, e_init), e_init])
    range = lepton_range(e_init, abs(injector.ν_pdg)==16)
    # Randomly sample zenith uniform in phase space
    θ = acos(rand(Uniform(cos(injector.θmax), cos(injector.θmin))))
    # Randomly sample azimuth
    ϕ = rand(Uniform(injector.ϕmin, injector.ϕmax))
    d = Direction(θ, ϕ)
    # Sample impact parameter uniformly on a disc
    b = injector.r_injection .* sqrt(rand())
    # Sample angle on disc 
    ψ = rand(Uniform(0, 2π))
    # Rotate to plane perpendicular to event direction
    p_near = SVector{3}(perpendicular_plane(θ, ϕ, b, ψ))
    p_int = sample_interaction_vertex(
        injector,
        p_near,
        d,
        range,
        geo
    )
    initial_state = Particle(
        injector.ν_pdg,
        e_init,
        p_int,
        d,
        nothing
    )
    final_state = Particle(
        injector.ν_pdg  - sign(injector.ν_pdg),
        e_final,
        p_int,
        d,
        initial_state
    )
    # Reverse Direction since Track tells us where we're going
    # But Particle direction tells us where it is from
    texit = Track(p_int, reverse(d), geo.box)
    ranges = computeranges(texit, geo)
    proposal_result = propagate(final_state, ranges)
    # Pass PROPOSAL output to CORSIKA
    event = InjectionEvent(initial_state, final_state, proposal_result)
    return event
end