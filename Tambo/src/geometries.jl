defaults = (zup=5units.km, zdown=50units.km, depth1=12units.km, depth2=21.4units.km)

struct Coord{T<:Float64}
    latitude::T
    longitude::T
end

function Base.show(io::IO, coord::Coord)
    println(
        io,
        """
        latitude (degrees): $(coord.latitude * 180 / π)°
        longitude (degrees): $(coord.longitude * 180 / π)°
        """
    )
end

"""
    latlong_to_xy(lat, long, latmin, longmin)

Function to calculate the xy coordinate for a latitude and longitude
coordinate. This implicitly assumes that Δθ and Δϕ are small so that
the sphere is locally floating. All angles are in radians !!!!!
latmin and longmin are -15.73975004° and -72.336236836° for the current spline

# Example
```julia-repl
julia> latmin, longmin = deg2rad(-15.73975004), deg2rad(-72.336236836);

julia> tapay_lat, tapay_long = deg2rad-15.63), deg2rad(-72.16;

julia> x_tapay, y_tapay = latlong_to_xy(tapay_lat, tapay_long, latmin, longmin) ./ units.km
(18.861842060081564, 12.203647647037522)
```
"""
function latlong_to_xy(lat, long, latmin, longmin)
    r = 6_371.0 * units.km
    x = r * cos(longmin) * (lat - latmin)
    y = r * (long - longmin)
    return x, y
end

struct Direction
    θ::Float64
    ϕ::Float64
    proj::SVector{3,Float64}
end

function Direction(θ::T, ϕ::U) where {T,U<:Number}
    θ, ϕ = Float64(θ), Float64(ϕ)
    proj = SVector{3}([cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)])
    return Direction(θ, ϕ, proj)
end

function Direction(x::T, y::U, z::V) where {T,U,V<:Number}
    x, y, z = promote(x, y, z)
    proj = SVector{3}([x, y, z]) ./ norm((x, y, z))
    θ = acos(proj.z)
    ϕ = atan(proj.y, proj.x)
    return Direction(θ, ϕ, proj)
end

function Direction(sv::SVector{3})
    return Direction(sv.x, sv.y, sv.z)
end

function Direction(ip::SVector{3}, fp::SVector{3})
    return Direction(ip .- fp)
end

function Direction(pp_vector::PyObject)
    return Direction(pp_vector.x, pp_vector.y, pp_vector.z)
end

function Base.show(io::IO, d::Direction)
    print(
        io,
        """
        θ (degrees): $(d.θ * 180 / π)°
        ϕ (degrees): $(d.θ * 180 / π)°
        proj: [$(proj.x), $(proj.y), $(proj.z)]
        """
    )
end

Base.reverse(d::Direction) = Direction(-d.proj...)

Base.:/(sv::SVector{3}, d::Direction) = sv ./ d.proj
Base.:*(m, d::Direction) = m .* d.proj

struct Box{T<:Number}
    c1::SVector{3,T}
    c2::SVector{3,T}
end

function Box(c1, c2)
    c1, c2 = promote(c1, c2)
    c1, c2 = SVector{3}(c1), SVector{3}(c2)
    return Box(c1, c2)
end

function Box(c)
    c1 = [0, 0, 0]
    c2 = c
    return Box(c1, c2)
end

function Box(x, y, z)
    c = [x, y, z]
    return Box(c)
end

function inside(sv::SVector{3}, b::Box)
    e1 = b.c1
    e2 = b.c2
    s = sign.(e2 .- e1)
    is_in = all(s .* e1 .< s .* sv < s .* e2)
    return is_in
end

inside(x, y, z, f::Function) = z < f(x, y)
inside(sv::SVector{3}, f::Function) = inside(sv.x, sv.y, sv[3], f)

function load_spline(p, key="spline")
    f = jldopen(p, "r")
    spl = f[key]
    return spl
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

function Geometry(spl::Dierckx.Spline2D, xyzoffset::SVector{3}, depths::Vector, ρs::Vector)
    zup = defaults.zup
    zdown = defaults.zdown
    knots = spl.tx, spl.ty
    xmin, xmax = minimum(knots[1]) * units[:m], maximum(knots[1]) * units[:m]
    ymin, ymax = minimum(knots[2]) * units[:m], maximum(knots[2]) * units[:m]
    xyzmin = SVector{3}([xmin - xyzoffset.x, ymin - xyzoffset.y, -zdown])
    xyzmax = SVector{3}([xmax - xyzoffset.x, ymax - xyzoffset.y, zup])
    valley(x, y) = valley_helper(x, y, xyzoffset, spl)
    b = Box(xyzmin, xyzmax)
    zboundaries = -xyzoffset.z .- depths
    return Geometry(valley, b, xyzoffset, units.ρair0, units.ρrock0, zboundaries, ρs)
end

function Geometry(spl::Dierckx.Spline2D, depths::Vector, ρs::Vector)
    knots = spl.tx, spl.ty
    xmin, xmax = minimum(knots[1]) * units[:m], maximum(knots[1]) * units[:m]
    ymin, ymax = minimum(knots[2]) * units[:m], maximum(knots[2]) * units[:m]
    xmid, ymid = (xmax - xmin) / 2, (ymax - ymin) / 2
    xyzoffset = SVector{3}([
        xmid, ymid, evaluate(spl, xmid / units[:m], ymid / units[:m]) * units[:m]
    ])
    return Geometry(spl, xyzoffset, depths, ρs)
end

function Geometry(spl::Dierckx.Spline2D, xyoffset::SVector{2})
    depths = []
    ρs = []
    return Geometry(spl, xyoffset, depths, ρs)
end

function Geometry(spl::Dierckx.Spline2D)
    depths = []
    ρs = []
    return Geometry(spl, depths, ρs)
end

function Geometry(spl::Dierckx.Spline2D, xyoffset::SVector{2}, depths::Vector, ρs::Vector)
    z = spl(xyoffset.x / units[:m], xyoffset.y / units[:m]) * units[:m]
    xyzoffset = SVector{3}([xyoffset.x, xyoffset.y, z])
    return Geometry(spl, xyzoffset, depths, ρs)
end

function Geometry(spl_path::String)
    spl = load_spline(spl_path)
    return Geometry(spl)
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
    return evaluate(valley_spl, xm, ym) * units[:m] - xyzoffset[3]
end

function (g::Geometry)(x, y)
    return g.valley(x, y)
end
