using LinearAlgebra:norm
using StaticArrays
using JLD2
using Dierckx

include("units.jl")

defaults = (
    zup=5units.km,
    zdown=50units.km,
    depth1=12units.km,
    depth2=21.4units.km
)

struct Coord{T<:Float64}
    lat::T
    long::T
end

function get_distance_coordinate(latitude, longitude, latmin, longmin)
    latmid = (latitude + latmin)/2.0
    m_per_deg_lon = (111132.954 - (559.822 * cos( 2.0 * latmid )) + (1.175 * cos( 4.0 * latmid)) + (0.0023 * cos( 6.0 * latmid)))
    m_per_deg_lat = (111412.82 * cos(latmid)) - (93.5*cos(latmid*3)) + (0.118*cos(5*latmid))
    delta_lat = latitude - latmin 
    delta_long = longitude - longmin 
    x = delta_long * (m_per_deg_lon * 180/pi)
    y = delta_lat * (m_per_deg_lat * 180/pi)
    x, y
end

struct Direction
    θ::Float64
    ϕ::Float64
    proj::SVector{3, Float64}
end

function Direction(θ::T, ϕ::U) where {T,U<:Number}
    θ, ϕ = Float64(θ), Float64(ϕ)
    proj = SVector{3}([cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)])
    Direction(θ, ϕ, proj)
end

function Direction(x::T, y::U, z::V) where {T,U,V<:Number}
    x, y, z = promote(x, y, z)
    proj = SVector{3}([x,y,z]) ./ norm((x,y,z))
    θ = acos(proj.z)
    ϕ = atan(proj.y, proj.x)
    Direction(θ, ϕ, proj)
end

function Direction(sv::SVector{3})
    Direction(sv.x, sv.y, sv.z)
end

function Direction(ip::SVector{3}, fp::SVector{3})
    Direction(ip .- fp)
end

Base.:/(sv::SVector{3}, d::Direction) = sv ./ d.proj
Base.:*(m, d::Direction) = m.*d.proj

struct Box{T<:Number}
    c1::SVector{3, T}
    c2::SVector{3, T}
end

function Box(c1, c2)
    c1, c2 = promote(c1, c2)
    c1, c2 = SVector{3}(c1), SVector{3}(c2)
    Box(c1, c2)
end

function Box(c)
    c1 = [0,0,0]
    c2 = c
    Box(c1, c2)
end

function Box(x, y, z)
    c = [x,y,z]
    Box(c)
end

function inside(sv::SVector{3}, b::Box)
    e1 = b.c1
    e2 = b.c2
    s = sign.(e2 .- e1)
    is_in = all(s .* e1 .< s .* sv < s .* e2)
    is_in
end

inside(x, y, z, f::Function) = z < f(x, y)
inside(sv::SVector{3}, f::Function) = inside(sv.x, sv.y, sv[3], f)

function load_spline(p, key="spline")
    f = jldopen(p, "r")
    spl = f[key]
    spl
end

struct Geometry
    valley::Function
    box::Box
    tambo_offset::SVector{3}
    ρair::Float64
    ρrock::Float64
    zboundaries::Vector{Float64}    
    ρs::Vector{Float64}
end

function Geometry(
    spl::Dierckx.Spline2D,
    xyzoffset::SVector{3},
    depths::Vector,
    ρs::Vector,
    # zdown=defaults.zdown,
    # zup=defaults.zup
)
    zup = defaults.zup
    zdown = defaults.zdown
    knots = spl.tx, spl.ty
    xmin, xmax = minimum(knots[1]) * units[:m], maximum(knots[1]) * units[:m]
    ymin, ymax = minimum(knots[2]) * units[:m], maximum(knots[2]) * units[:m]
    xyzmin = SVector{3}([xmin-xyzoffset.x, ymin-xyzoffset.y, -zdown])
    xyzmax = SVector{3}([xmax-xyzoffset.x, ymax-xyzoffset.y, zup])
    valley(x, y) = valley_helper(x, y, xyzoffset, spl)
    b = Box(xyzmin, xyzmax)
    zboundaries = -xyzoffset.z .- depths
    Geometry(valley, b, xyzoffset, units.ρair0, units.ρrock0, zboundaries, ρs)
end

function Geometry(
    spl::Dierckx.Spline2D,
    depths::Vector,
    ρs::Vector,
    # zdown=defaults.zdown,
    # zup=defaults.zup,
)
    knots = spl.tx, spl.ty
    xmin, xmax = minimum(knots[1]) * units[:m], maximum(knots[1]) * units[:m]
    ymin, ymax = minimum(knots[2]) * units[:m], maximum(knots[2]) * units[:m]
    xmid, ymid = (xmax - xmin)/2, (ymax - ymin)/2
    xyzoffset = SVector{3}([
        xmid,
        ymid,
        evaluate(spl, xmid / units[:m], ymid / units[:m]) * units[:m]
    ])
    Geometry(spl, xyzoffset, depths, ρs)
end

function Geometry(
    spl::Dierckx.Spline2D,
    xyoffset::SVector{2},
    # zdown=defaults.zdown,
    # zup=defaults.zup, 
)
    depths = []
    ρs = []
    Geometry(spl, xyoffset, depths, ρs)
end

function Geometry(
    spl::Dierckx.Spline2D,
    # zdown=defaults.zdown,
    # zup=defaults.zup,
)
    depths = []
    ρs = []
    Geometry(spl, depths, ρs)
end

function Geometry(
    spl::Dierckx.Spline2D,
    xyoffset::SVector{2},
    depths::Vector,
    ρs::Vector,
    # zdown=defaults.zup,
    # zup=defaults.zdown, 
)
    z = spl(xyoffset.x / units[:m], xyoffset.y / units[:m]) * units[:m]
    xyzoffset = SVector{3}([xyoffset.x, xyoffset.y, z])
    Geometry(spl, xyzoffset, depths, ρs)
end

function Geometry(spl_path::String)
    spl = load_spline(spl_path)
    Geometry(spl)
end

"""
    valley_helper(x, y, tc, valley_spl)

Function for calling the Python `valley_spl`. This converts `x` and `y` from the
TAMBO-centered coordinate system which uses natural units to the spline
coordinate system which uses meters. It then converts the spline result back to
natural units. The offset between these coordinate systems is given by `tc`, the
center of TAMBO in the spline coordinate system.

# Example
```julia-repl
julia> v = valley_helper(
           1000m, 4000m, 
           TPoint(1.08e11, 9.255e10, 1.35e10), 
           spl
       ) / m

       488.4400118544429
```
"""
function valley_helper(x, y, xyzoffset, valley_spl)
    # Translate to spline coordinate system
    xt, yt = x + xyzoffset[1], y + xyzoffset[2]
    # Convert to meters
    xm, ym = xt / units[:m], yt / units[:m]
    # Evaluate spline and convert back to natural units
    evaluate(valley_spl, xm, ym)*units[:m] - xyzoffset[3]
end

function (g::Geometry)(x, y)
    g.valley(x,y)
end
