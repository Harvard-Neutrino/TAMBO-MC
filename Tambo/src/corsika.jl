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

function should_do_corsika(event::ProposalResult, geo::Geometry)
    # Check if going right direction
    norm_pos = event.propped_state.position ./ norm(event.propped_state.position)
    # TODO I made up this Number. We should check if it makes any sense
    if event.propped_state.direction.proj'norm_pos < -0.1
        return false
    end
    return has_unobstructed_path(event.propped_state.position, geo)
#=
function should_do_corsika(event::ProposalResult, geo::Geometry)
    # Check if going right direction
    decay_pos = event.decay_products[1].position/units.m
    
    """
    If the decay occurs inside the mountain, cut. 
    """
    if inside(decay_pos,geo) 
    return false 
    end 

    propped_dir = event.propped_state.direction
    
    distance, point, dot = intersect(dstate.position,state.direction,plane) 
    
    """
    The particle has to travel backwards to reach the plane. 
    We don't want particles that have to travel backwards to reach the plane. 
    """
    if distance/units.m < 0
    return false 
    end 
    
    """
    The particle direction is in the same direction as the plane's normal vector. 
    To intercept with positive distance, the particle would have to be traveling from inside the mountain.  
    """

    if dot > 0 
    return false 
    end 

    """
    Cutting near-orthogonal particle directions with the plane normal. 
    """
    if abs(dot) < 1e-3
    return false  
    end

    """
    point[3] = z-intercept of particle and TAMBO plane.If the elevation in TAMBO coords is greater than 10km, cut. 
    """
    if point[3]/units.m > 10000
    return false 
    end 

    """
    If the distance length is greater than 20km between particle position and intercept with TAMBO plane, cut. 
    """

    if distance/units.m > 20000
    return false 
    end 

    return has_unobstructed_path(event.decay_state[1].position, geo)
=#

end

