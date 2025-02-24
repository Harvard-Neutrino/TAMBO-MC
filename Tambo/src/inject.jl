struct Injector
    nu_pdg::Int
    powerlaw::PowerLaw
    xs::CrossSection
    anglesampler::UniformAngularSampler
    injectionshape::AbstractInjectionShape
    geo::Geometry
end

struct InjectionEvent
    event_id::Int
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

    if "r_injection" in keys(config) && "l_endcap" in keys(config)
        injectionshape = SymmetricInjectionCylinder(
            config["r_injection"] * units.m,
            config["l_endcap"] * units.m
        )
    elseif "length" in keys(config) && "width" in keys(config)
        injectionshape = InjectionPlane(
            config["length"] * units.m,
            config["width"] * units.m,
            geo.tambo_normal.proj
        )
    else
        error("Unknown options for loading injection shape")
    end
    return Injector(config["nu_pdg"], pl, xs, anglesampler, injectionshape, geo)
end

function inject_event(injector::Injector, event_id::Int, tr_seed::Int)
    event = inject_event(
        event_id,
        injector.nu_pdg,
        injector.powerlaw,
        injector.xs,
        injector.anglesampler,
        injector.injectionshape,
        injector.geo,
        tr_seed
    )
    return event
end

function sample_interaction_vertex(
    endcap_distance::Float64,
    closest_approach::SVector{3},
    d::Direction,
    range::Float64,
    geo::Geometry
)
    track = Track(closest_approach, reverse(d), geo.box)
    segments = computesegments(track, geo)
    #@show segments
    tot_X = endcapcolumndepth(track, endcap_distance, range, segments)
    #@show tot_X
    X = rand(Uniform(0.0, tot_X))
    λ_int = inversecolumndepth(track, X, geo, segments)
    p_int = track(λ_int)
    return p_int, tot_X
end

#function sample_interaction_vertex(
#    plane::InjectionPlane,
#    plane_position::SVector{3},
#    d::Direction,
#    range::Float64,
#    geo::Geometry
#)
#    track = Track(plane_position, reverse(d), geo.box)
#    segments = computesegments(track, geo)
#    tot_X = minimum([range, totalcolumndepth(track, segments)])
#    #if tot_X==0
#    #    @show range
#    #    @show segments
#    #    @show plane_position
#    #end
#    X = rand(Uniform(0.0, tot_X))
#    λ_int = inversecolumndepth(track, X, geo, segments)
#    p_int = track(λ_int)
#    return p_int, tot_X
#end

"""
    endcapcolumndepth(t::Track, l_endcap::Float64, range::Float64, segments::Vector{Segment})

TBW
"""
function endcapcolumndepth(
    t::Track,
    l_endcap::Float64,
    range::Float64,
    segments::Vector{Segment}
)
    cd = totalcolumndepth(t, segments)
    #@show cd
    if t.norm <= l_endcap
        return cd
    end
    cd_endcap = minimum([columndepth(t, l_endcap / t.norm, segments) + range, cd])
    #@show l_endcap / t.norm
    #@show columndepth(t, l_endcap / t.norm, segments)
    #@show cd_endcap
    return cd_endcap
end

function determine_injection_start(
    closest_approach::SVector{3},
    boundary::SVector{3},
    geo::Geometry
)
    track = Track(closest_approach, boundary)
    segments = computesegments(track, geo)
    if length(segments)==1
        return closest_approach, false
    end

    start = closest_approach
    saw_air, no_endcap = false, false
    for segment in segments
        if segment.medium_name=="Air"
            if ~saw_air
                start = segment.pstart
            end
            saw_air = true
        elseif saw_air && segment.medium_name=="StandardRock"
            start = segment.pstart
            l_endcap = 0.0
            break
        end
    end
    return start, no_endcap
end

"""

TBW
"""
function inject_event(
    event_id::Int,
    ν_pdg::Int,
    power_law::PowerLaw,
    xs::CrossSection,
    anglesampler::UniformAngularSampler,
    cylinder::SymmetricInjectionCylinder,
    geo::Geometry,
    tr_seed::Int
)
    # Randomly sample direction
    direction = Direction(rand(anglesampler)...)
    # Rotation to plane perpindicular to direction
    rotator = (RotX(direction.θ) * RotZ(π / 2 - direction.ϕ))'
    # Find point of closest approach
    closest_approach = rotator * rand(cylinder)
    # Find where the particle entered the TAMBO region
    xb = intersect(closest_approach, reverse(direction), geo.box)

    # Sample initial neutrino energy
    proposed_e_init = rand(power_law)
    # Make particle and do tau regenerations
    proposed_particle = Particle(ν_pdg, proposed_e_init, xb, direction, nothing)
    particle_entry, physX = tr_propagate(proposed_particle, geo.tambo_offset.z, tr_seed)
    # Set energy of neutrino when enters the box
    e_final = rand(xs, particle_entry.energy)

    # Find where to start counting CD from
    #injection_start, no_endcap = determine_injection_start(closest_approach, xb, geo)
    injection_start, no_endcap = closest_approach, false
    l_endcap = cylinder.l_endcap
    if no_endcap
        l_endcap = 0.0
    end
    range = lepton_range(e_final, ν_pdg - sign(ν_pdg))
    if range > 0
        p_int, genX = sample_interaction_vertex(l_endcap, closest_approach, direction, range, geo)
    else
        p_int, genX = closest_approach, floatmin()
    end
    final_state = Particle(ν_pdg - sign(ν_pdg), e_final, p_int, direction, particle_entry)
    # Now our generation column depth is the same as our physical column depth so we set them equal
    event = InjectionEvent(
        event_id,
        particle_entry,
        proposed_particle,
        final_state,
        genX,
        genX
    )
    return event
end

#function inject_event(
#    ν_pdg::Int,
#    power_law::PowerLaw,
#    xs::CrossSection,
#    anglesampler::UniformAngularSampler,
#    plane::InjectionPlane,
#    geo::Geometry,
#    tr_seed::Int
#)
#    direction = Direction(rand(anglesampler)...)
#    # Rotation to plane perpindicular to direction
#    #rotator = (RotX(direction.θ) * RotZ(π / 2 - direction.ϕ))'
#    #closest_approach = rotator * rand(injectionshape)
#    plane_position = rand(plane)
#    xb = intersect(plane_position, reverse(direction), geo.box)
#
#    proposed_e_init = rand(power_law)
#    proposed_particle = Particle(ν_pdg, proposed_e_init, xb, direction, nothing)
#    particle_entry, physX = tr_propagate(proposed_particle, geo.tambo_offset.z, tr_seed)
#    e_final = rand(xs, particle_entry.energy)
#
#    range = lepton_range(e_final, ν_pdg - sign(ν_pdg))
#    if range==0
#        @show particle_entry
#        p_int, genX = plane_position, 0.0
#    else
#        p_int, genX = sample_interaction_vertex(plane, plane_position, direction, range, geo)
#    end
#    final_state = Particle(ν_pdg - sign(ν_pdg), e_final, p_int, direction, particle_entry)
#    # Now our generation column depth is the same as our physical column depth so we set them equal
#    event = InjectionEvent(
#        particle_entry,
#        proposed_particle,
#        final_state,
#        genX,
#        genX
#    )
#    return event
#end
