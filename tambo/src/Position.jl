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