"""
CORSIKA orgin is at intersection of particle and observation plane
x-axis is aligned north
y-axis is aligned west
Look at page 97 of the CORISKA manual for this information
"""

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

#struct CorsikaInfo
#    run_number::Int
#    particle::Particle
#    observation_location::SVector{3}
#    plane::Plane
#    map::AffineMap
#end
#
#function CorsikaInfo(particle::Particle, geo::Geometry, run_number::Int, plane::Plane)
#
#    t = Track(particle.position, particle.direction, geo.box)
#    obs_loc = intersect(t, plane) .- particle.position
#    mapping = tambo_to_corsika_mapping(particle, geo)
#
#    return CorsikaInfo(
#        run_number,
#        particle,
#        obs_loc,
#        plane,
#        mapping
#    )
#end
#
#function Base.show(io::IO, info::CorsikaInfo)
#    print(
#        io, 
#        """
#        run_number: $(info.run_number)
#        energy: $(info.energy / units.GeV) GeV
#        direction: $(info.direction)
#        observation_height: $(info.observation_height / units.m) m
#        plane: $(info.plane)
#        initial_elevation: $(info.initial_elevation / units.m) m"""
#    )
#end
#
#"""
#    tambo_to_corsika(d::Direction)
#
#TBW
#"""
#function tambo_to_corsika(d::Direction)
#    return Direction(π - d.θ, d.ϕ - π / 2)
#end
#
#"""
#    corsika_to_tambo(d::Direction)
#
#
#"""
#function corsika_to_tambo(d::Direction)
#    return Direction(π - d.θ, d.ϕ + π / 2)
#end
#
#function tambo_to_corsika_mapping(particle::Particle, geo::Geometry)
#    """
#    The corsika coordinate system has its origin at the x- and y-positions
#    of the initial particle, and its z-position at sea-level
#    """
#    corsika_origin = SVector{3}([
#        particle.position.x,
#        particle.position.y,
#        -geo.tambo_offset.z
#    ])
#    # The CORSIKA coordinate system aligns the x-axis with magnetic North
#    # and the y-axis with West. We align 
#    rot = RotZ(π / 2)
#    trans = Translation(corsika_origin)
#    mapping = inv(trans ∘ rot)
#    return mapping
#end
