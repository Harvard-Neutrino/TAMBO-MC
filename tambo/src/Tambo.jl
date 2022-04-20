module Tambo

push!(LOAD_PATH, @__DIR__)

using Geometry: GenerationRegion, TPoint, sample
using Tracks: Track, Direction, reverse, total_column_depth
using Particles: Particle
using PowerLaws: PowerLaw, sample

mutable struct TAMBOSim
    n::Int
    gr::GenerationRegion
    ν_pdg::Int
    γ::Float64
    emin::Float64
    emax::Float64
    θmin::Float64
    θmax::Float64
    ϕmin::Float64
    ϕmax::Float64
    function TAMBOSim()
        n = 0
        gr = GenerationRegion("/Users/jlazar/research/TAMBO-MC/resources/tambo_spline.npy")
        ν_pdg = 16
        γ = 2
        emin = 1e6
        emax = 1e9
        θmin = 0
        θmax = π
        ϕmin = 0
        ϕmax = 2π
        new(n, gr, ν_pdg, γ, emin, emax, θmin, θmax, ϕmin, ϕmax)
    end
end

function (ts::TAMBOSim)()
    # sample initial neutrino energies according to PL
    pl = PowerLaw(ts.γ, ts.emin, ts.emax)
    νs = [Particle(ts.ν_pdg, pl) for _ in 1:ts.n]
    ## Sample τ energies according to physical distribution
    #τ_pdg = ν_pdg+1*sign(ν_pdg)
    #τs = Particle.(Ref(t_pdg), tau_energy.(νs))
    ts, cds = make_tracks(ts)
end

function tau_energy(eν)
    u = rand()
    decay_param(u, eν)
end

function tau_energy(ν::Particle)
    u = rand()
    decay_param(u, ν.energy)
end

function make_tracks(
    n::Int,
    gr::GenerationRegion,
    θmin::Float64,
    θmax::Float64,
    ϕmin::Float64,
    ϕmax::Float64
)
    xyz = sample(n, gr)
    # Sample zenith angles equally in cos(theta)
    # TOOD actually use the args
    θ = acos.(rand(n).*2.0.-1)
    # Sample azimuths
    ϕ = rand(n) .* (ϕmax-ϕmin) .+ ϕmin
    # This track starts in the middle of the box and then exits
    ts = Track.(TPoint.(xyz), Direction.(θ, ϕ), Ref(gr.box))
    # We reverse it to get a particle that is entering
    ts = reverse.(ts)
    cd = total_column_depth.(ts, Ref(gr))
    #=
    We obviously need some sort of event structure for this, but I'll leave it
    for now because I want to merge this stuff and get on with my day.....
    =#
    ts, cd
end

function make_tracks(ts::TAMBOSim)
    make_tracks(ts.n, ts.gr, ts.θmin, ts.θmax, ts.ϕmin, ts.ϕmax)
end

function sample_properties(
    n::Int,
    gr::GenerationRegion,
    ν_pdg::Int,
    γ::Float64,
    emin::Float64,
    emax::Float64,
    θmin::Float64,
    θmax::Float64,
    ϕmin::Float64,
    ϕmax::Float64
)
    eτ = tau_energy.(eν)
    τ = Particle.(Ref(ν_pdg+1*sign(ν_pdg)), eτ)
    #=
    Sample points where the τ decays uniformly within the box. 
    We could also make a different box for each particle 
    and optimize the size to make the simulation faster a la LeptonInjector
    =#
    xyz = sample(n, gr)
    # Sample zenith angles equally in cos(theta)
    θ = acos.(rand(n).*2.0.-1)
    # Sample azimuths
    ϕ = rand(n) .* (ϕmax-ϕmin) .+ ϕmin
    # This track starts in the middle of the box and then exits
    t = Track.(TPoint.(xyz), Direction.(θ, ϕ), Ref(gr.box))
    # We reverse it to get a particle that is entering
    t = reverse.(t)
    cd = total_column_depth.(t, Ref(gr.valley))
    #= 
    We obviously need some sort of event structure for this, but I'll leave it
    for now because I want to merge this stuff and get on with my day.....
    =#
    ν, τ, t, cd
end

end # module
