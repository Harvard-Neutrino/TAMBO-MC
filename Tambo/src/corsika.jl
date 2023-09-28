struct CorsikaEvent
  pdg::Int
  kinetic_energy::Number
  x::Float64
  y::Float64
  z::Float64
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
  return CorsikaEvent(pdg, kinetic_energy, x, y, z, t, weight)
end

struct DirectionMap
    ϕ_offset::Float64
end

function (dm::DirectionMap)(direction::Direction)
    return Direction(π - direction.θ, direction.ϕ + dm.ϕ_offset)
end

function Base.inv(dm::DirectionMap)
    return DirectionMap(-dm.ϕ_offset)
end

struct PositionMap
    # I think that the origin should go, but w.e. for now
    corsika_origin::SVector{3}
    mapping::AffineMap
end

function PositionMap(particle::Particle, geo::Geometry)
    # The corsika coordinate system has its origin at the x- and y-positions
    # of the initial particle, and its z-position at sea-level
    corsika_origin = SVector{3}([
        particle.position.x,
        particle.position.y,
        -geo.tambo_offset.z
    ])
    # The CORSIKA coordinate system aligns the x-axis with magnetic North
    # and the y-axis with West. We align 
    rot = LinearMap(RotZ(π / 2))
    trans = Translation(corsika_origin)
    mapping = inv(trans ∘ rot)
    return PositionMap(corsika_origin, mapping)
end

function (pm::PositionMap)(x)
    return pm.mapping(x)
end

function Base.inv(pm::PositionMap)
    return PositionMap(pm.corsika_origin, inv(pm.mapping))
end

struct CorsikaMap
    direction_map::DirectionMap
    position_map::PositionMap
end

function CorsikaMap(particle::Particle, geo::Geometry)
    pm = PositionMap(particle, geo)
    dm = DirectionMap(-π / 2)
    return CorsikaMap(dm, pm)
end

function (cm::CorsikaMap)(x::SVector{3})
    return cm.position_map(x)
end

function (cm::CorsikaMap)(d::Direction)
    return cm.direction_map(d)
end

function Base.inv(cm::CorsikaMap)
    return CorsikaMap(inv(cm.direction_map), inv(cm.position_map))
end

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
    point[3] = z-intercept of particle and TAMBO plane.If the elevation in TAMBO coords is greater than 10km, cut. 
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

function check_track_length(event, plane, geo; verbose=false)
    """
    If the distance length is greater than 20km between particle position and intercept with TAMBO plane, cut. 
    """
    decay_pos = event.propped_state.position
    propped_dir = event.propped_state.direction
    #plane = Tambo.Plane(minesite_normal_vec, minesite_coord, geo)      
    distance, _, _ = intersect(decay_pos, propped_dir, plane) 
    if distance > 20000 * units.km
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
    #plane = Tambo.Plane(minesite_normal_vec, minesite_coord, geo)      
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

#=
function should_do_corsika(event::ProposalResult, geo::Geometry)
    # Check if going right direction
    norm_pos = event.propped_state.position ./ norm(event.propped_state.position)
    # TODO I made up this Number. We should check if it makes any sense
    if event.propped_state.direction.proj'norm_pos < -0.1
        return false
    end
    return has_unobstructed_path(event.propped_state.position, geo)
=#
function should_do_corsika(event::ProposalResult, plane::Plane, geo::Geometry, criteria::Array{Function}; verbose=false, check_mode=false)
    # Seems like we should do `check_mode` via dependency injection but weeeeeeeeeeee
    if check_mode
        b = Bool[]
    end

    for criterion in criteria
        v = criterion(event,plane, geo; verbose=verbose)
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
        check_intersections
    ]
    return should_do_corsika(event, plane, geo, checks; verbose=verbose, check_mode=check_mode)
end
