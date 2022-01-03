
function get_distance_coordinate(latitude, longitude, latmin, longmin)
    latmid = (latitude + latmin)/2.0
    m_per_deg_lon = (111132.954 - (559.822 * cos( 2.0 * latmid )) + (1.175 * cos( 4.0 * latmid)) + (0.0023 * cos( 6.0 * latmid)))
    m_per_deg_lat = (111412.82 * cos(latmid)) - (93.5*cos(latmid*3)) + (0.118*cos(5*latmid))
    delta_lat = latitude - latmin 
    delta_long = longitude - longmin 

    x = delta_long * (m_per_deg_lon * 180/pi)
    y = delta_lat * (m_per_deg_lat * 180/pi)

    x,y
end

function construct_point(latitude, longitude, elevation; latmin = 0.0, longmin = 0.0)
    x,y = get_distance_coordinate(latitude,longitude,latmin, longmin)
    Point(longitude, latitude, x, y, elevation)
end

struct Point{T<:Real}
    longitude::T
    latitude::T
    x::T
    y::T
    z::T
end

function construct_point(phi, theta)
    x = cos(theta)*sin(phi)
    y = sin(theta)*sin(phi)
    z = cos(phi)
    Direction(phi,theta,x,y,z)
end

struct Direction{T<:Real}
    phi::T
    theta::T
    x::T
    y::T
    z::T
end