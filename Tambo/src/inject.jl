struct Injector
    nu_pdg::Int
    powerlaw::PowerLaw
    xs::CrossSection
    anglesampler::UniformAngularSampler
    injectionvolume::SymmetricInjectionCylinder
    geo::Geometry
end

struct InjectionEvent
    entry_state::Particle
    initial_state::Particle
    final_state::Particle
    physX::Float64
    genX::Float64
end

function Injector(config::Dict, geo::Geometry)
    pl = PowerLaw(
        config["gamma"],
        config["emin"] * units.GeV,
        config["emax"] * units.GeV
    )
    xs = CrossSection(
        config["xs_dir"],
        config["xs_model"],
        config["nu_pdg"],
        config["interaction"]
    )
    anglesampler = UniformAngularSampler(
        deg2rad(config["thetamin"]),
        deg2rad(config["thetamax"]),
        deg2rad(config["phimin"]),
        deg2rad(config["phimax"])
    )

    injectionvolume = SymmetricInjectionCylinder(
        config["r_injection"] * units.m,
        config["l_endcap"] * units.m
    )
    return Injector(config["nu_pdg"], pl, xs, anglesampler, injectionvolume, geo)
end

function inject_event(injector::Injector, tr_seed::Int)
    event = inject_event(
        injector.nu_pdg,
        injector.powerlaw,
        injector.xs,
        injector.anglesampler,
        injector.injectionvolume,
        injector.geo,
        tr_seed
    )
    return event
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

TBW
"""
function inject_event(
    ν_pdg::Int,
    power_law::PowerLaw,
    xs::CrossSection,
    anglesampler::UniformAngularSampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry,
    tr_seed::Int
)
    direction = Direction(rand(anglesampler)...)
    # Rotation to plane perpindicular to direction
    rotator = (RotX(direction.θ) * RotZ(π / 2 - direction.ϕ))'
    closest_approach = rotator * rand(injectionvolume)
    xb = intersect(closest_approach, reverse(direction), geo.box)

    proposed_e_init = rand(power_law)
    proposed_particle = Particle(ν_pdg, proposed_e_init, xb, direction, nothing)
    
    particle_entry, physX = tr_propagate(proposed_particle, geo.tambo_offset.z, tr_seed)
    
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
