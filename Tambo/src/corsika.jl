Base.@kwdef mutable struct CORSIKAConfig
    parallelize_corsika::Bool = false
    thinning::Float64 = 1e-6 
    hadron_ecut::Float64 = 0.05units.GeV
    em_ecut::Float64 = 0.01units.GeV
    photon_ecut::Float64 = 0.002units.GeV
    mu_ecut::Float64 = 0.05units.GeV 
    shower_dir::String = "showers"
    singularity_path::String = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-env.simg"
    corsika_path::String = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/corsika"
    corsika_sbatch_path::String = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/TAMBO-MC/scripts/corsika_parallel.sbatch"

    tambo_coordinates::Coord = whitepaper_coord
    plane_orientation::Direction = whitepaper_normal_vec
    proposal_events::Vector{ProposalResult} = [] 
end

struct CORSIKAPropagator
    config::CORSIKAConfig 
    geo::Geometry 
end 

function (propagator::CORSIKAPropagator)(; track_progress=true)
    n = length(propagator.config.proposal_events)
    proposal_events = propagator.config.proposal_events
    iter = 1:n 
    if track_progress 
        iter = ProgressBar(iter) 
    end 

    geo = propagator.geo 
    #plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)
    plane = Tambo.Plane(propagator.config.plane_orientation, propagator.config.tambo_coordinates, geo)
    indices = []
    for (proposal_idx,proposal_event) in enumerate(proposal_events)
        update(iter)
        if should_do_corsika(proposal_event,plane,geo)
            for (decay_idx,decay_event) in enumerate(proposal_event.decay_products)

                #checks if decay product is a neutrino 
                #skips if is 
                #wanted to keep indices lined up so checking one at at ime
                if check_neutrino(decay_event)
                    push!(indices,[proposal_idx,decay_idx])
                else 
                    continue 
                end 

                if propagator.config.parallelize_corsika 
                    continue 
                else 
                    corsika_run(decay_event,propagator,proposal_idx,decay_idx; parallelize_corsika=false)
                end 
            end
        end
    end
    
    if propagator.config.parallelize_corsika 
        if track_progress 
            n = length(indices)
            println("Running CORSIKA in parallel for $n showers")
        end 
        corsika_parallel(proposal_events,propagator,indices)
    end 
    return indices 
end 
   


# function corsika_inject_event(injector::Injector)
#     return event
# end

# function (injector::Injector)(; track_progress=true)
#     seed!(injector.config.seed)
#     iter = 1:injector.config.n
#     if track_progress
#         iter = ProgressBar(iter)
#     end
#     return [corsika_inject_event(injector) for _ in iter]
# end

function corsika_parallel(proposal_events,propagator,indices)

    if length(indices) > 10000
        println("Slow down there partner...")
        println("FASRC forbids more than 10k concurrent jobs!!")
    end 

    for i in indices
        
        proposal_idx = i[1]
        decay_idx = i[2]
        
        corsika_run(proposal_events[proposal_idx].decay_products[decay_idx],
        propagator, 
        proposal_idx,
        decay_idx; 
        parallelize_corsika=true)
        
    end 
    #parallelize_corsika_exec = `sbatch $sbatch_dir`
    #run(corsika_exec);
end 

function corsika_run(
    pdg::Int64,
    energy::Float64,
    zenith::Float64, 
    azimuth::Float64, 
    inject_pos::SVector{3},
    intercept_pos::SVector{3}, 
    plane::SVector{3}, 
    obs_z::Float64, 
    thinning::Float64, 
    ecuts::SVector{4},
    singularity_path::String,
    corsika_path::String,
    corsika_sbatch_path::String,
    outdir::String,
    proposal_index::Int64,
    decay_index::Int64; 
    parallelize_corsika=parallelize_corsika
    )
    rawinject_x,rawinject_y,rawinject_z = inject_pos
    x_intercept,y_intercept,z_intercept = intercept_pos 
    xdir,ydir,zdir = plane 
    
    #convert to CORSIKA internal units of GeV
    emcut,photoncut,mucut,hadcut = ecuts/units.GeV 
    total_index = string(proposal_index) *"_"* string(decay_index)

    if parallelize_corsika 
        corsika_parallel_exec = "singularity exec $singularity_path $corsika_path --pdg $pdg --energy $energy --zenith $zenith --azimuth $azimuth --xpos $rawinject_x --ypos $rawinject_y --zpos $rawinject_z -f $outdir/shower_$total_index --xdir $xdir --ydir $ydir --zdir $zdir --observation-height $obs_z --force-interaction --x-intercept $x_intercept --y-intercept $y_intercept --z-intercept $z_intercept --emcut $emcut --photoncut $photoncut --mucut $mucut --hadcut $hadcut --emthin $thinning"
        run(`sbatch $corsika_sbatch_path $corsika_parallel_exec`)
    else 
        corsika_exec = `singularity exec $singularity_path $corsika_path --pdg $pdg --energy $energy --zenith $zenith --azimuth $azimuth --xpos $rawinject_x --ypos $rawinject_y --zpos $rawinject_z -f $outdir/shower_$total_index --xdir $xdir --ydir $ydir --zdir $zdir --observation-height $obs_z --force-interaction --x-intercept $x_intercept --y-intercept $y_intercept --z-intercept $z_intercept --emcut $emcut --photoncut $photoncut --mucut $mucut --hadcut $hadcut --emthin $thinning`
        run(corsika_exec)
    end 
end 

function corsika_run(decay_event::Particle,propagator::CORSIKAPropagator,proposal_idx::Int64,decay_idx::Int64; parallelize_corsika=parallelize_corsika)
    geo = propagator.geo
    thinning = propagator.config.thinning
    ecuts = SVector{4}([propagator.config.em_ecut,propagator.config.photon_ecut,propagator.config.mu_ecut,propagator.config.hadron_ecut])
    outdir = propagator.config.shower_dir 
    singularity_path = propagator.config.singularity_path
    corsika_path = propagator.config.corsika_path
    corsika_sbatch_path = propagator.config.corsika_sbatch_path
    return corsika_run(decay_event::Particle,geo,thinning,ecuts,singularity_path,corsika_path,corsika_sbatch_path,outdir,proposal_idx::Int64,decay_idx::Int64; parallelize_corsika=parallelize_corsika)
end

function corsika_run(decay_event::Particle,geo::Geometry,thinning::Float64,ecuts,singularity_path::String,corsika_path::String,corsika_sbatch_path::String,outdir::String,proposal_idx::Int64,decay_idx::Int64; parallelize_corsika=parallelize_corsika)
    plane = Plane(whitepaper_normal_vec, whitepaper_coord, geo)

    pdg = decay_event.pdg_mc
    energy = decay_event.energy/units.GeV
    zenith = decay_event.direction.θ
    azimuth = decay_event.direction.ϕ
    corsika_map = CorsikaMap(decay_event,geo)
    obs_z = corsika_map(plane.x0)[3] / units.km

    distance, point, dot = intersect(decay_event.position, decay_event.direction, plane)
    distance = distance/units.km
    intercept_pos = point/units.km 
    inject_pos = decay_event.position/units.km 
  
    return corsika_run(
        pdg,
        energy,
        zenith,
        azimuth,
        inject_pos,
        intercept_pos,
        plane.n̂.proj,
        obs_z, 
        thinning,
        ecuts,
        singularity_path,
        corsika_path,
        corsika_sbatch_path,
        outdir,
        proposal_idx,
        decay_idx;
        parallelize_corsika=parallelize_corsika
        )
end 

struct CorsikaEvent
  pdg::Int
  kinetic_energy::Number
  pos::SVector{3}
  time::Float64
  weight::Float64
end

function CorsikaEvent(a::PyObject, particle_idx::Int)
  pdg = a["pdg"].to_numpy()[particle_idx]
  t = a["time"].to_numpy()[particle_idx]
  x = a["x"].to_numpy()[particle_idx]
  y = a["y"].to_numpy()[particle_idx]
  z = a["z"].to_numpy()[particle_idx]
  kinetic_energy = a["kinetic_energy"].to_numpy()[particle_idx]
  weight = a["weight"].to_numpy()[particle_idx]
  return CorsikaEvent(pdg, kinetic_energy, SVector{3}([x, y, z]), t, weight)
end

function CorsikaEvent(weight)
    return CorsikaEvent(0, 0.0, 0.0, 0.0, 0.0, 0.0, weight)
end
#
# This is not really true... Fix this to include events and not just indices
struct Hit
  mod::DetectionModule
  event::CorsikaEvent
end

#struct DirectionMap
#    ϕ_offset::Float64
#end
#
#function (dm::DirectionMap)(direction::Direction)
#    return Direction(π - direction.θ, direction.ϕ + dm.ϕ_offset)
#end
#
#function Base.inv(dm::DirectionMap)
#    return DirectionMap(-dm.ϕ_offset)
#end
#
#struct PositionMap
#    # I think that the origin should go, but w.e. for now
#    corsika_origin::SVector{3}
#    mapping::AffineMap
#end
#
#function PositionMap(particle::Particle, geo::Geometry)
#    # The corsika coordinate system has its origin at the x- and y-positions
#    # of the initial particle, and its z-position at sea-level
#    corsika_origin = SVector{3}([
#        particle.position.x,
#        particle.position.y,
#        -geo.tambo_offset.z
#    ])
#    # The CORSIKA coordinate system aligns the x-axis with magnetic North
#    # and the y-axis with West. We align 
#    rot = LinearMap(RotZ(π / 2))
#    trans = Translation(corsika_origin)
#    mapping = inv(trans ∘ rot)
#    return PositionMap(corsika_origin, mapping)
#end
#
#function (pm::PositionMap)(x)
#    return pm.mapping(x)
#end
#
#function Base.inv(pm::PositionMap)
#    return PositionMap(pm.corsika_origin, inv(pm.mapping))
#end
#
#struct CorsikaMap
#    direction_map::DirectionMap
#    position_map::PositionMap
#end
#
#function CorsikaMap(particle::Particle, geo::Geometry)
#    pm = PositionMap(particle, geo)
#    dm = DirectionMap(-π / 2)
#    return CorsikaMap(dm, pm)
#end
#
#function (cm::CorsikaMap)(x::SVector{3})
#    return cm.position_map(x)
#end
#
#function (cm::CorsikaMap)(d::Direction)
#    return cm.direction_map(d)
#end
#
#function Base.inv(cm::CorsikaMap)
#    return CorsikaMap(inv(cm.direction_map), inv(cm.position_map))
#end

function check_inside_mtn(event, plane, geo; verbose=false)
    decay_pos = event.propped_state.position
    if inside(decay_pos, geo) 
        if verbose
            println("inside mountain")
        end
        return false 
    end 
    return true 
end

function check_right_direction(event, plane, geo; verbose=false)
    """
    The particle has to travel backwards to reach the plane. 
    We don't want particles that have to travel backwards to reach the plane. 
    """
    decay_pos = event.propped_state.position
    propped_dir = event.propped_state.direction
    #plane = Tambo.Plane(minesite_normal_vec, minesite_coord, geo)      
    distance, _, _ = intersect(decay_pos, propped_dir, plane) 
    if distance/units.m < 0
        if verbose
            println("negative distance")
        end
        return false 
    end
    return true 
end

function check_plane_dot(event, plane, geo; verbose=false)
    """
    The particle direction is in the same direction as the plane's normal vector. 
    To intercept with positive distance, the particle would have to be traveling from inside the mountain.  
    """

    decay_pos = event.propped_state.position
    propped_dir = event.propped_state.direction
    #plane = Tambo.Plane(minesite_normal_vec, minesite_coord, geo)      
    _, _, dot = intersect(decay_pos, propped_dir, plane) 
    if dot > 0 
        if verbose
            println("DOT > 0")
        end
        return false 
    end 
    return true 
end

function check_near_orthogonal(event, plane, geo; verbose=false)
    """
    Cutting near-orthogonal particle directions with the plane normal. 
    """
    decay_pos = event.propped_state.position
    propped_dir = event.propped_state.direction
    #plane = Tambo.Plane(minesite_normal_vec, minesite_coord, geo)      
    _, _, dot = intersect(decay_pos, propped_dir, plane) 
    if abs(dot) < 1e-3
        if verbose
            println("ORTHOGONAL")
        end
        return false  
    end
    return true 
end

function check_z_intercept(event, plane, geo; verbose=false)
    """
    point[3] = z-intercept of particle and TAMBO plane. If the elevation in TAMBO coords is greater than 10km, cut. 
    """
    decay_pos = event.propped_state.position
    propped_dir = event.propped_state.direction
    #plane = Tambo.Plane(minesite_normal_vec, minesite_coord, geo)      
    _, point, _ = intersect(decay_pos, propped_dir, plane) 
    if point.z > 10 * units.km
        if verbose
            println("Z-intercept GREATER THAN 10km")
        end
        return false 
    end 
    return true 
end

function check_track_length(event, plane, geo; verbose=false, max_distance=20units.km)
    """
    If the distance length is greater than 20km between particle position and intercept with TAMBO plane, cut. 
    """
    decay_pos = event.propped_state.position
    propped_dir = event.propped_state.direction
    #plane = Tambo.Plane(minesite_normal_vec, minesite_coord, geo)      
    distance, _, _ = intersect(decay_pos, propped_dir, plane) 
    if distance > max_distance
        if verbose
            println("Distance length greater than 20km") 
        end
        return false 
    end 
    return true 
end

function check_intersections(event, plane, geo; verbose=false)
    decay_pos = event.propped_state.position
    propped_dir = event.propped_state.direction
    _, point, _ = intersect(decay_pos, propped_dir, plane) 
    t = Track(decay_pos, point)
    intersections = intersect(t, geo)

    """
    If the intersection with the mountain > 1 then we cut.
    When it intersects TAMBO plane.
    """
    if length(intersections) > 1
        return false
    end
    
    return true 
end

function check_passed_through_rock(event::ProposalResult, plane, geo; verbose=false, thresh=4units.km)
    track = Tambo.Track(
        event.continuous_losses.position,
        reverse(event.propped_state.direction),
        geo.box
    )
    segments = Tambo.computesegments(track, geo)
    
    ℓ_inrock = 0
    for segment in segments
        if segment.medium_name=="Air"
            continue
        end
        ℓ_inrock += segment.length
    end
    return ℓ_inrock >= thresh
end

function check_neutrino(daughter_particle::Particle; verbose = false)
    pdg = daughter_particle.pdg_mc
    if abs(pdg) == 16 || abs(pdg) == 14 || abs(pdg) == 12 
        if verbose 
            println("daughter particle is a neutrino")
        end 
        return false
    else 
        if verbose 
            println("daughter is not a neutrino")
        end
        return true 
    end 
end

function should_do_corsika(event::ProposalResult, plane::Plane, geo::Geometry, criteria::Array{Function}; verbose=false, check_mode=false)
    # Seems like we should do `check_mode` via dependency injection but weeeeeeeeeeee
    if check_mode
        b = Bool[]
    end

    for criterion in criteria
        v = criterion(event, plane, geo; verbose=verbose)
        if check_mode
            push!(b, v)
        else
            if !v
                return false
            end
        end
    end

    if check_mode
        return b
    else
        return true
    end
end

function should_do_corsika(event::ProposalResult, plane::Plane, geo::Geometry; verbose=false, check_mode=false)
    checks = [
        check_inside_mtn,
        check_right_direction,
        check_plane_dot,
        check_near_orthogonal,
        check_z_intercept,
        check_track_length,
        check_intersections,
        check_passed_through_rock,
    ]
    return should_do_corsika(event, plane, geo, checks; verbose=verbose, check_mode=check_mode)
end
