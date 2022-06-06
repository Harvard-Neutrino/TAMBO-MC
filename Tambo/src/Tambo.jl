module Tambo

export TAMBOSim

include("powerLaws.jl")
include("tracks.jl")
include("units.jl")

using StaticArrays
using Rotations
using Random

mutable struct TAMBOSim
    n::Int
    geo::Geometry
    ν_pdg::Int
    γ::Float64
    #emin::Quantity{Float64, Unitful.𝐋^2*Unitful.𝐌 /Unitful.𝐓^2}
    #emax::Quantity{Float64, Unitful.𝐋^2*Unitful.𝐌 /Unitful.𝐓^2}
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
  
    function TAMBOSim()
        n = 0
        geo = Geometry("/Users/jlazar/research/TAMBO-MC/resources/tambo_spline.jld2")
        ν_pdg = 16
        γ = 2
        emin = 1e6units[:GeV]
        emax = 1e9units[:GeV]
        pl = nothing
        θmin = 0
        θmax = π
        ϕmin = 0
        ϕmax = 2π
        r_injection = 900units[:m]
        l_endcap = 1units[:km]
        seed = 0
        new(n, geo, ν_pdg, γ, emin, emax, pl, θmin, θmax, ϕmin, ϕmax, r_injection, l_endcap, seed)
    end
end

function (ts::TAMBOSim)()
    Random.seed!(ts.seed)
    verify_ts!(ts)
    inject_events(ts)
end

function verify_ts!(ts::TAMBOSim)
    change_pl = false
    if change_pl || ==(ts.pl, nothing)
        ts.pl = PowerLaw(ts.γ, ts.emin, ts.emax)
    end
    if !(abs(ts.ν_pdg) ∈ [14, 16])
        throw(ErrorException, "invalid ν_pdg. must be in [±14, ±16]")
    end
end

"""
    muon_range(e)

Compute the 99.9% column depth for a μ with energy `e` using parametrization
from https://doi.org/10.1016/j.cpc.2021.108018

# Example
```julia-repl
julia> muon_range(1 PeV)
9.837812593932223e7 kg m⁻²
```
"""
function muon_range(e)
    da = 0.1777/units[:GeV]/units[:mwe]
    db = 2.09*10^-4/units[:mwe]
    cd = log(1 + e * da / db) / db
end

"""
    tau_range(e)

Compute the 99.9% column depth a secondary μ from primary τ decay with energy `e`
using parametrization from https://doi.org/10.1016/j.cpc.2021.108018

# Example
```julia-repl
julia> tau_range(1 PeV)
9.905162093356338e7 kg m⁻²
```
"""
function tau_range(e)
    da = 4.7e-13/units[:GeV]/units[:mwe]
    db = 2.63e-5/units[:mwe]
    cd = log(1 + e * da / db) / db + muon_range(e)
end

"""
    lepton_range(e, ν_pdg)

Copute the 99.9% column_depth of a charged lepton emerging from a ν CC event
using the parametrization from https://doi.org/10.1016/j.cpc.2021.108018

# Example
```julia-repl
julia> lepton_range(1 PeV, 16)
9.905162093356338e7 kg m⁻²

julia> lepton_range(1 PeV, 14)
9.837812593932223e7 kg m⁻²
```
"""
function lepton_range(e, ν_pdg)
    if abs(ν_pdg==16)
        range = tau_range(e)
    elseif abs(ν_pdg==14)
        range = muon_range(e)
    end
    range
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

#function sample_column_depth(ti::Track, to::Track, ts::TAMBOSim, range)
#    cdi = total_column_depth(ti, ts.geo.valley)
#    cdo = total_column_depth(to, ts.geo.valley)
#    if ti.norm <= ts.l_endcap
#        #println("If you're seeing this a lot the injection region is too big")
#        cdi_endcap = cdi
#    else
#        cdi_endcap = minimum(
#            [column_depth(ti, ts.l_endcap/ti.norm, ts.geo.valley) + range, cdi]
#        )
#    end
#    if to.norm <= ts.l_endcap
#        #println("If you're seeing this a lot the injection region is too big")
#        cdo_endcap = cdo
#    else
#        cdo_endcap = column_depth(to, ts.l_endcap/to.norm, ts.geo.valley)
#    end
#    cd = rand() * (cdo_endcap + cdi_endcap)
#    cd < cdi_endcap ? tr = ti : tr = to
#    cd = abs(cdi_endcap - cd)
#    cd, tr
#end

function sample_column_depth(t::Track, ts::TAMBOSim, range)
    tot_X = total_column_depth(t, ts.geo.valley)
    if ti.norm <= ts.l_endcap
        cdi_endcap = cdi
    else
        cdi_endcap = minimum(
            [column_depth(ti, ts.l_endcap/ti.norm, ts.geo.valley) + range, cdi]
        )
    end
    if to.norm <= ts.l_endcap
        cdo_endcap = cdo
    else
        cdo_endcap = column_depth(to, ts.l_endcap/to.norm, ts.geo.valley)
    end
    cd = rand() * (cdo_endcap + cdi_endcap)
    cd < cdi_endcap ? tr = ti : tr = to
    cd = abs(cdi_endcap - cd)
    cd, tr
end

function sample_tau_energy(eν, νtype, xs)

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

function inject_events(ts::TAMBOSim)
    
    # Sample an energy
    e = rand(ts.n, ts.pl)
    range = lepton_range.(e, Ref(ts.ν_pdg))
    # Randomly sample zenith uniform in phase space
    θ = acos.(rand(ts.n) .* (cos(ts.θmin)-cos(ts.θmax)) .+ cos(ts.θmax))
    # Randomly sample azimuth
    ϕ = rand(ts.n) .* (ts.ϕmax-ts.ϕmin) .+ ts.ϕmin
    # Sample impact parameter uniformly on a disc
    b = ts.r_injection .* sqrt.(rand(ts.n))
    # Sample angle on disc 
    ψ = rand(ts.n) .* 2π
    # Rotate to plane perpendicular to event direction
    p_near = SVector{3}.(perpendicular_plane.(θ, ϕ, b, ψ))
    # Make track from point of closest approach to point of entry
    ti = Track.(p_near, Direction.(θ, ϕ), Ref(ts.geo.box))
    # Make track from point of closest approach to point of exit
    to = Track.(p_near, Direction.(π.-θ, mod.(ϕ.+π, 2π)), Ref(ts.geo.box))
    ipoint = intersect.(ti, Ref(ts.geo.box))
    fpoint = intersect.(to, Ref(ts.geo.box))
    tr = Track.(ipoint, fpoint)
    ixs = intersect.(tr, ts.geo.valley)
    tot_cd = total_column_depth.(tr, Ref(ts.geo.valley), ixs)
    cd = rand(ts.n) .* tot_cd
    # Find affine parameter where we have traversed proper column depth
    λ_int = inverse_column_depth.(tr, cd, ts.geo.valley, ixs)
    # Convert affine parameter to a physical location
    # TODO This feels wrong but I don't know what is right
    p_int = [tr[i](λ_int[i]) for i in eachindex(tr)]
    # Sample an outgoing lepton energy

    ## Make a list of media that the lepton sample_properties
    ## Pass to Jorge's function
    ## Pass PROPOSAL output to CORSIKA
    Event.(e, θ, ϕ, b, ψ, ti, to, tot_cd, p_int, p_near, tr, λ_int)
end

end # module
