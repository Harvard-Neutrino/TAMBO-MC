module Position

export Coord, 
       TPoint,
       Direction

#-------------

struct Coord{T<:Float64}
    lat::T
    long::T
end

#-------------

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

using LinearAlgebra
struct TPoint
    x
    y
    z
    function TPoint(c::Coord, elevation; lat_min = 0.0, long_min = 0.0)
        x,y = get_distance_coordinate(c.lat, c.long,lat_min, long_min)
        new(x, y, elevation)
    end
    function TPoint(x, y, z)
        new(x, y, z)
    end
    function TPoint(phi, theta)
        x = cos(phi)*sin(theta)
        y = sin(phi)*sin(theta)
        z = cos(theta)
        new(x, y, z)
    end
end


Base.:+(p1::TPoint, p2::TPoint) = TPoint(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
Base.:-(p1::TPoint, p2::TPoint) = TPoint(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
Base.:/(p1::TPoint, p2::TPoint) = TPoint(p1.x/p2.x, p1.y/p2.y, p1.z/p2.z)
Base.:/(p1::TPoint, t)          = TPoint(p1.x/t, p1.y/t, p1.z/t)
Base.:*(p1::TPoint, p2::TPoint) = TPoint(p1.x*p2.x, p1.y*p2.y, p1.z*p2.z)
Base.:*(t, p::TPoint)           = TPoint(t*p.x, t*p.y, t*p.z)
Base.:*(p::TPoint, t)           = t*p
LinearAlgebra.norm(p::TPoint)    = norm((p.x, p.y, p.z))
#-------------

struct Direction
    phi
    theta
    point::TPoint
    function Direction(phi, theta, x, y, z)
        len   = norm((x,y,z))
        x     = x/len
        y     = y/len
        z     = z/len
        new(phi, theta, TPoint(x, y, z))
    end
    function Direction(x, y, z)
        len   = norm((x,y,z))
        x     = x/len
        y     = y/len
        z     = z/len
        phi   = atan(y,x)
        theta = acos(z)
        new(phi, theta, TPoint(x, y, z))
    end
    function Direction(phi, theta)
        x = cos(theta)*sin(phi)
        y = sin(theta)*sin(phi)
        z = cos(phi)
        new(phi, theta, TPoint(x, y, z))
    end
end

Base.:+(d::Direction, p::TPoint)  = d.point+p
Base.:+(p::TPoint, d::Direction)  = d+p
Base.:-(d::Direction, p::TPoint)  = d.point-p
Base.:-(p::TPoint, d::Direction)  = p-d.point
Base.:*(d::Direction, t::Float64) = d.point*t
Base.:*(t::Float64, d::Direction) = d*t
Base.:*(d::Direction, p::TPoint)  = d.point*p
Base.:*(p::TPoint, d::Direction)  = d*p
Base.:/(d::Direction, p::TPoint)  = d.point/p
Base.:/(p::TPoint, d::Direction)  = p/d.point

end
