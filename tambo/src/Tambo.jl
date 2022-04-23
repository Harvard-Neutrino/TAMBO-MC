module Tambo

push!(LOAD_PATH, @__DIR__)

using Geometry: GenerationRegion, TPoint, sample
using Tracks
using Particles: Particle
using PowerLaws
using Units
using StaticArrays
using Rotations
using Unitful

mutable struct TAMBOSim
    n::Int
    gr::GenerationRegion
    ν_pdg::Int
    γ::Float64
    emin::Quantity{Float64, Unitful.𝐋^2*Unitful.𝐌 /Unitful.𝐓^2}
    emax::Quantity{Float64, Unitful.𝐋^2*Unitful.𝐌 /Unitful.𝐓^2}
    pl::Union{PowerLaw, Nothing}
    θmin::Float64
    θmax::Float64
    ϕmin::Float64
    ϕmax::Float64
    r_injection::Quantity{Float64, Unitful.𝐋}
    l_endcap::Quantity{Float64, Unitful.𝐋}

    function TAMBOSim()
        n = 0
        gr = GenerationRegion("/Users/jlazar/research/TAMBO-MC/resources/tambo_spline.npy")
        ν_pdg = 16
        γ = 2
        emin = 1e6GeV
        emax = 1e9GeV
        pl = nothing
        θmin = 0
        θmax = π
        ϕmin = 0
        ϕmax = 2π
        r_injection = 900m
        l_endcap = 1km
        new(n, gr, ν_pdg, γ, emin, emax, pl, θmin, θmax, ϕmin, ϕmax, r_injection, l_endcap)
    end
end

function verify_ts!(ts::TAMBOSim)
    change_pl = false
    #if !=(ts.pl, nothing)
    #    if ts.emin != ts.pl.emin
    #        @test_logs (:warn, "Power law lower limits don't match. Will redefine")
    #        change_pl = true
    #    end
    #    if ts.emax != ts.pl.emax
    #        change_pl = true
    #        @test_logs (:warn, "Power law upper limits don't match. Will redefine")
    #    end
    #end
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
    da = 0.1777/GeV/mwe
    db = 2.09*10^-4/mwe
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
    println(e)
    da = 4.7e-13/GeV/mwe
    db = 2.63e-5/mwe
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
function perpendicular_plane(θ, ϕ, b, ψ; return_transform=false)
    # Construct vector in the plane of normal coordinate system
    bv = SVector{3}([b*cos(ψ), b*sin(ψ), b*0])
    # Make matrix to rotate to perpendicular plane
    r = (Rotations.RotX(θ) * RotZ(π/2-ϕ))'
    # Rotate to perpendicular plane
    if return_transform
        return r * bv, r
    else
        return r * bv
    end
end

function sample_column_depth(ti::Track, to::Track, ts::TAMBOSim, range)
    cdi = total_column_depth(ti, ts.gr.valley)
    cdo = total_column_depth(to, ts.gr.valley)
    if ti.norm <= ts.l_endcap
        cdi_endcap = cdi
    else
        cdi_endcap = minimum(
            [column_depth(ti, ts.l_endcap/ti.norm, ts.gr.valley) + range, cdi]
        )
    end
    if to.norm <= ts.l_endcap
        cdo_endcap = cdo
    else
        cdo_endcap = column_depth(to, ts.l_endcap/to.norm, ts.gr.valley)
    end
    cd = rand() * (cdo_endcap + cdi_endcap)
    cd < cdi_endcap ? tr = ti : tr = to
    cd = abs(cdi_endcap - cd)
    cd, tr
end


function inject_events(ts::TAMBOSim)
    # Sample an energy
    e = rand(ts.pl)
    range = lepton_range(e, ts.ν_pdg)
    # Randomly sample zenith uniform in phase space
    θ = acos(rand() * (cos(ts.θmin)-cos(ts.θmax)) +cos(ts.θmax))
    # Randomly sample azimuth
    ϕ = rand() * (ts.ϕmax-ts.ϕmin) + ts.ϕmin
    # Sample impact parameter uniformly on a disc
    b = ts.r_injection * sqrt.(rand())
    # Sample angle on disc 
    ψ = rand() * 2π
    # Rotate to plane perpendicular to event direction
    # This is the point of closest approach
    p_near = TPoint(perpendicular_plane(θ, ϕ, b, ψ))
    # Make track from point of closest approach to point of entry
    ti = Track(p_near, Direction(θ, ϕ), ts.gr.box)
    # Make track from point of closest approach to point of exit
    to = Track(p_near, Direction(π-θ, mod(ϕ+π, 2\pi)), ts.gr.box)
    # Calculate the total column seen on the way in and way out
    cd, tr = sample_column_depth(ti, to, ts, range)
    # Find affine parameter where we have traversed proper column depth
    λ_int = inverse_column_depth(tr, cd, ts.gr.valley)
    # Convert affine parameter to a physical location
    p_int = tr(λ_int)
    # Sample an outgoing lepton energy
    # Make a list of media that the lepton sample_properties
    # Pass to Jorge's function
    # Pass PROPOSAL output to CORSIKA
    Event(revelant_stuff)
end

end # module
#function inject_events(ts::TAMBOSim)
#    # Randomly sample azimuth
#    ϕ = rand(ts.n) .* (ts.ϕmax-ts.ϕmin) .+ ts.ϕmin
#    # Randomly sample zenith uniform in phase space
#    θ = acos.(rand(ts.n) .* (cos(ts.θmin)-cos(ts.θmax))) .+ ts.θmin
#    # Sample impact parameter uniformly on a disc
#    b = ts.r_injection .* sqrt.(rand(ts.n))
#    # Sample angle on disc 
#    ψ = rand(ts.n) .* 2π
#    # Rotate to plan perpendicular to 
#end

#
#function (ts::TAMBOSim)()
#    # sample initial neutrino energies according to PL
#    pl = PowerLaw(ts.γ, ts.emin, ts.emax)
#    νs = [Particle(ts.ν_pdg, pl) for _ in 1:ts.n]
#    ## Sample τ energies according to physical distribution
#    #τ_pdg = ν_pdg+1*sign(ν_pdg)
#    #τs = Particle.(Ref(t_pdg), tau_energy.(νs))
#    ts, cds = make_tracks(ts)
#end
#
#function tau_energy(eν)
#    u = rand()
#    decay_param(u, eν)
#end
#
#function tau_energy(ν::Particle)
#    u = rand()
#    decay_param(u, ν.energy)
#end
#
#function make_tracks(
#    n::Int,
#    gr::GenerationRegion,
#    θmin::Float64,
#    θmax::Float64,
#    ϕmin::Float64,
#    ϕmax::Float64
#)
#    xyz = sample(n, gr)
#    # Sample zenith angles equally in cos(theta)
#    # TOOD actually use the args
#    θ = acos.(rand(n).*2.0.-1)
#    # Sample azimuths
#    ϕ = rand(n) .* (ϕmax-ϕmin) .+ ϕmin
#    # This track starts in the middle of the box and then exits
#    ts = Track.(TPoint.(xyz), Direction.(θ, ϕ), Ref(gr.box))
#    # We reverse it to get a particle that is entering
#    ts = reverse.(ts)
#    cd = total_column_depth.(ts, Ref(gr))
#    #=
#    We obviously need some sort of event structure for this, but I'll leave it
#    for now because I want to merge this stuff and get on with my day.....
#    =#
#    ts, cd
#end
#
#function make_tracks(ts::TAMBOSim)
#    make_tracks(ts.n, ts.gr, ts.θmin, ts.θmax, ts.ϕmin, ts.ϕmax)
#end
#
#function sample_properties(
#    n::Int,
#    gr::GenerationRegion,
#    ν_pdg::Int,
#    γ::Float64,
#    emin::Float64,
#    emax::Float64,
#    θmin::Float64,
#    θmax::Float64,
#    ϕmin::Float64,
#    ϕmax::Float64
#)
#    eτ = tau_energy.(eν)
#    τ = Particle.(Ref(ν_pdg+1*sign(ν_pdg)), eτ)
#    #=
#    Sample points where the τ decays uniformly within the box. 
#    We could also make a different box for each particle 
#    and optimize the size to make the simulation faster a la LeptonInjector
#    =#
#    xyz = sample(n, gr)
#    # Sample zenith angles equally in cos(theta)
#    θ = acos.(rand(n).*2.0.-1)
#    # Sample azimuths
#    ϕ = rand(n) .* (ϕmax-ϕmin) .+ ϕmin
#    # This track starts in the middle of the box and then exits
#    t = Track.(TPoint.(xyz), Direction.(θ, ϕ), Ref(gr.box))
#    # We reverse it to get a particle that is entering
#    t = reverse.(t)
#    cd = total_column_depth.(t, Ref(gr.valley))
#    #= 
#    We obviously need some sort of event structure for this, but I'll leave it
#    for now because I want to merge this stuff and get on with my day.....
#    =#
#    ν, τ, t, cd
#end