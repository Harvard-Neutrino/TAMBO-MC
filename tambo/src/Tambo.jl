# This is mostly sudocode at this point. Don't take it too seriously
module tambo

push!(LOAD_PATH, @__DIR__)

using Geometry: GenerationRegion, TPoint
using Tracks: Track, Direction, reverse, total_column_depth
using Particles: Particle

geo = load_geo("geo.json")
gr = GenerationRegion(spline_path)

function sample_properties(n, ν_pdg)
    #=
    We need to sample the direction first because that dictates the relationship
    between thei energies of the primary and secondary neutrinos
    =#
    # We only want to sample over upgoing events 
    θ = acos.(rand(n).-1)
    # Sample over all azimuths
    ϕ = rand(n).*(2*π)
    # Sample primary energies according to physical distribution
    eν_primary = sample_power_law(n, gamma)
    # This is where the spline that Ibrahim is making comes in
    eν_final = secondary_energy.(eν_primary, θ)
    # Sample τ energy from decay distribution
    eτ = tau_energy.(eν_final)

    ν = Particle(ν_pdg, eν)
    τ = Particle(ν_pdg+1*sign(ν_pdg), eτ)
    #=
    Sample points where the τ decays uniformly within the box. 
    We could also make a different box for each particle 
    and optimize the size to make the simulation faster a la LeptonInjector
    =#
    xyz = [rand(3).*abs.(gr.box.c1.-gr.box.c2).+gr.box.minxyz for _ in 1:n]
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
