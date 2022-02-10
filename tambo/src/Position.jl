module Position

using LinearAlgebra
using Rotations
using CoordinateTransformations
using StaticArrays

export Coord, 
       TPoint,
       Box,
       CoordinateTransform

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
Base.:/(p1::TPoint, t) = TPoint(p1.x/t, p1.y/t, p1.z/t)
Base.:*(p1::TPoint, p2::TPoint) = TPoint(p1.x*p2.x, p1.y*p2.y, p1.z*p2.z)
Base.:*(t, p::TPoint) = TPoint(t*p.x, t*p.y, t*p.z)
Base.:*(p::TPoint, t) = t*p
LinearAlgebra.norm(p::TPoint) = norm((p.x, p.y, p.z))
Base.:*(r::Rotations.Rotation, p::TPoint) = TPoint(r*[p.x, p.y, p.z]...)
function (t::Translation)(p::TPoint)
    TPoint(t([p.x, p.y, p.z])...)
end


"""
Construct the `CoordinateTransform` object for transforming between coordinate systems.
These coordinate systems are related by a rotation and a translation. The new coordinate system
is first rotated about the z-axis (ψ), then the y-axis (θ), then the x-axis (ϕ).
After this it is translated by (x,y,z)
"""
struct CoordinateTransform
    translation::Translation
    rotation::RotXYZ
    function CoordinateTransform(x::SVector{3}, r::SVector{3})
        trans = Translation(x)
        rot = RotXYZ(r...)
        new(trans, rot)
    end
    function CoordinateTransform(x, r)
        x = SVector{3}(x)
        r = SVector{3}(r...)
        CoordinateTransform(x, r)
    end
    function CoordinateTransform(x, y, z, ϕ, θ, ψ)
        x = SVector{3}(x,y,z)
        r = SVector{3}(ϕ, θ, ψ)
        CoordinateTransform(x, r)
    end
end

function (ct::CoordinateTransform)(x)
    x = ct.rotation * x
    ct.translation(x)
end

function Base.inv(ct::CoordinateTransform)
    f(x) = inv(ct.rotation)*inv(ct.translation)(x)
end

struct Box
    # User configurable points
    lwh::SVector{3}
    transform::CoordinateTransform
    function Box(lwh)
        lwh = SVector{3}(lwh)
        # Make a null transform
        transform = CoordinateTransform(0,0,0,0,0,0)
        new(lwh, transform)
    end
    function Box(lwh, transform::CoordinateTransform)
        if any(lwh.<=0)
            error("The length width and height of the box must be strictly positive")
        end
        lwh = SVector{3}(lwh)
        new(lwh, transform)
    end
end

"Function to convert coordinates in a box to a common refernce frame"
universal_coordinates(x, b::Box) = b.transform(x)
"Function to convert coordinates in a common refernce frame to those in a box"
box_coordinates(x, b::Box) = inv(b.transform)(x)
function box_to_box(x, b1::Box, b2::Box)
    x = universal_coordinates(x, b1)
    x = box_coordinates(x, b2)
end

function is_inside(p::TPoint, box::Box)
    # Convert to the coordinate system of the box
    p = box_coordinates(p, box)
    is_in = true
    fieldnames = [:x, :y, :z]
    for idx in 1:length(fieldnames)
        fn = fieldnames[idx]
        # Get the lenght width or height of the box
        dim = box.lwh[idx]/2
        #=
        Center of the box is the origin of this coordinate system
        If you are more than half that away, then you are outside the box
        =#
        is_in = is_in && (-dim < getfield(p, fn) < dim)
        #=
        If it is out in one dimension, we don't need to keep checking
        This probably offers next to no speed up but it's good habit
        I guess
        =#
        if !is_in
            return false
        end
    end
    is_in
end

function is_inside(x, y, z, f::Function)
    fz = f(x, y)
    if z<fz
        is_in = true
    else
        is_in = false
    end
    is_in
end

function is_inside(p::TPoint, f::Function)
    is_inside(p.x, p.y, p.z, f)
end

end