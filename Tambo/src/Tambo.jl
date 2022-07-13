module Tambo

export TAMBOSim

include("powerLaws.jl")
include("tracks.jl")
include("units.jl")
include("crosssections.jl")
include("particles.jl")

using StaticArrays
using Rotations
using Random
using StatsBase: sample

mutable struct TAMBOSim
    n::Int
    geo::Geometry
    ν_pdg::Int
    γ::Float64
    emin::Float64
    emax::Float64
    pl::Union{PowerLaw, Nothing}
    θmin::Float64
    θmax::Float64
    ϕmin::Float64
    ϕmax::Float64
    r_injection::Float64
    l_endcap::Float64
    seed::Int64
    diff_xs::DifferentialXS
  
    function TAMBOSim()
        n = 0
        geo = Geometry("$(@__DIR__)/../../resources/tambo_spline.jld2")
        ν_pdg = 16
        γ = 2
        emin = 1e6units[:GeV]
        emax = 1e9units[:GeV]
        pl = PowerLaw(2, 1e6units[:GeV], 1e9units[:GeV])
        θmin = 0
        θmax = π
        ϕmin = 0
        ϕmax = 2π
        r_injection = 900units[:m]
        l_endcap = 1units[:km]
        seed = 0
        diff_xs = DifferentialXS(
            "$(@__DIR__)/../../resources/cross_sections/tables/csms_differential.h5",
            ν_pdg
        )
        new(n, geo, ν_pdg, γ, emin, emax, pl, θmin, θmax, ϕmin, ϕmax, r_injection, l_endcap, seed, diff_xs)
    end
end

function (ts::TAMBOSim)()
    Random.seed!(ts.seed)
    verify_ts!(ts)
    [inject_event(ts) for _ in 1:ts.n]
end

function verify_ts!(ts::TAMBOSim)
    ts.pl = PowerLaw(ts.γ, ts.emin, ts.emax)
    abs(ts.ν_pdg) ∈ [14, 16] || error("invalid ν_pdg. must be in [±14, ±16]")
end

"""
    perpendicular_plane(θ, ϕ, b, ψ, [return_transform])

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
    # Make matrix to rotate to perpendicular plane
    # TODO This seems inefficient
    r = (Rotations.RotX(θ) * RotZ(π/2-ϕ))'
    r * bv
end


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

struct Event
    e::Float64
    θ::Float64
    ϕ::Float64
    impact_parameter::Float64
    ψ::Float64
    incoming_track::Track
    outgoing_track::Track
    column_depth::Float64
    interaction_vertex::SVector{3}
    # These are here for debugging. Will go away eventually
    p_near::SVector{3}
    tr::Track
    λ_int::Float64
end

function inject_event(ts::TAMBOSim)
    # Sample an energy
    e = rand(ts.pl)
    range = lepton_range(e, abs(ts.ν_pdg)==16)
    # Randomly sample zenith uniform in phase space
    θ = acos(rand() * (cos(ts.θmin)-cos(ts.θmax)) + cos(ts.θmax))
    # Randomly sample azimuth
    ϕ = rand() * (ts.ϕmax-ts.ϕmin) + ts.ϕmin
    # Sample impact parameter uniformly on a disc
    b = ts.r_injection .* sqrt(rand())
    # Sample angle on disc 
    ψ = rand() * 2π
    # Rotate to plane perpendicular to event direction
    p_near = SVector{3}(perpendicular_plane(θ, ϕ, b, ψ))
    # Make track from point of closest approach to point of entry
    ti = Track(p_near, Direction(θ, ϕ), ts.geo.box)
    # Make track from point of closest approach to point of exit
    to = Track(p_near, Direction(π-θ, mod(ϕ+π, 2π)), ts.geo.box)
    # Compute the intersection of each track with the mountain
    rangesi = computeranges(ti, ts.geo)
    rangeso = computeranges(to, ts.geo)
    # Compute the colum depth for both incoming and outgoing portions
    # TODO figure out why the same tracks give different cds...
    # Specifically this happens when you choose θmin = θmax = π
    cdi = endcapcolumndepth(ti, ts.l_endcap, range, rangesi)
    cdo = endcapcolumndepth(to, ts.l_endcap, 0.0, rangeso)
    # sample column depth uniformly and subtract incoming column depth
    cd = cdi - rand() * (cdi+cdo)
    # If the remainder is positive, you need to be in incoming track, else outgoing
    cd > 0 ? tr = ti : tr = to
    cd > 0 ? ranges = rangesi : ranges = rangeso
    # Get rid of potential negative
    cd = abs(cd)
    # Find affine parameter where we have traversed proper column depth
    λ_int = inversecolumndepth(tr, cd, ts.geo, ranges)
    # Convert affine parameter to a physical location
    p_int = tr(λ_int)
    # Sample an outgoing lepton energy
    e_τ = sample(ts.diff_xs, e)
    # Outgoing track
    texit = Track(p_int, Direction(π-θ, mod(ϕ+π, 2π)), ts.geo.box)
    # Pass to Jorge's function
    # Pass PROPOSAL output to CORSIKA
    Event(e, θ, ϕ, b, ψ, ti, to, cd, p_int, p_near, tr, λ_int)
    #Event(e, θ, ϕ, b, ψ, ti, to, cd, p_int, p_near, tr, λ_int, e_τ), texit
end


end # module