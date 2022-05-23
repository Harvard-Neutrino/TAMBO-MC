module Geometries

using LinearAlgebra, StaticArrays
using PyCall
using Units: m

export Geometry,
       TPoint,
       Direction,
       Box,
       is_inside,
       sample_xyz


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
    function TPoint(θ, ϕ)
        x = cos(ϕ)*sin(θ)
        y = sin(ϕ)*sin(θ)
        z = cos(θ)
        new(x, y, z)
    end
    function TPoint(sv::SVector{3})
        new(sv...)
    end
    function TPoint(v::Vector)
        new(v...)
    end
end

Base.:+(p1::TPoint, p2::TPoint) = TPoint(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
Base.:+(p::TPoint, sv::SVector{3}) = TPoint(p.x+sv[1], p.y+sv[2], p.z+sv[3])
Base.:-(p1::TPoint, p2::TPoint) = TPoint(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
Base.:-(p::TPoint, sv::SVector) = TPoint(p.x-sv[1], p.y-sv[2], p.z-sv[3])
Base.:-(sv::SVector{3}, p::TPoint) = -1*(p-sv)
Base.:/(p1::TPoint, p2::TPoint) = TPoint(p1.x/p2.x, p1.y/p2.y, p1.z/p2.z)
Base.:/(p1::TPoint, t) = TPoint(p1.x/t, p1.y/t, p1.z/t)
Base.:*(p1::TPoint, p2::TPoint) = TPoint(p1.x*p2.x, p1.y*p2.y, p1.z*p2.z)
Base.:*(t, p::TPoint) = TPoint(t*p.x, t*p.y, t*p.z)
Base.:*(p::TPoint, t) = Base.:*(t, p)
LinearAlgebra.norm(p::TPoint) = norm((p.x, p.y, p.z))

struct Direction
    θ
    ϕ
    x_proj
    y_proj
    z_proj
end

function Direction(θ, ϕ)
    θ, ϕ = Float64(θ), Float64(ϕ)
    x = cos(ϕ)*sin(θ)
    y = sin(ϕ)*sin(θ)
    z = cos(θ)
    Direction(θ, ϕ, x, y, z)
end

function Direction(x, y, z)
    x,y,z = (x,y,z)./norm((x,y,z))
    θ = acos(z)
    ϕ = atan(y, x)
    Direction(θ, ϕ, x, y, z)
end

function Direction(p::TPoint)
    Direction(p.x, p.y, p.z)
end

function Direction(ip::TPoint, fp::TPoint)
    Direction(ip-fp)
end

Base.:/(p::TPoint, d::Direction) = TPoint(p.x/d.x_proj, p.y/d.y_proj, p.z/d.z_proj)
Base.:*(m, d::Direction) = TPoint(m*d.x_proj, m*d.y_proj, m*d.z_proj)

struct Box
    c1::SVector{3}
    c2::SVector{3}
end

function Box(c1, c2)
    c1 = SVector{3}(c1)
    c2 = SVector{3}(c2)
    Box(c1, c2)
end

function Box(x, y, z)
    c1 = SVector{3}([0,0,0])
    c2 = SVector{3}([x,y,z])
    Box(c1, c2)
end

function Box(c)
    c1 = SVector{3}([0,0,0])
    c2 = SVector{3}(c)
    Box(c1, c2)
end

function sample(n::Int, b::Box)
    minxyz = SVector{3}(minimum.(zip(b.c1, b.c2)))
    [rand(3).*abs.(b.c1.-b.c2).+b.minxyz for _ in 1:n]
end

function sample(b::Box)
    minxyz = SVector{3}(minimum.(zip(b.c1, b.c2)))
    rand(3) .* abs.(b.c1 .- b.c2) .+ minxyz
end

function is_inside(p::TPoint, b::Box)
    is_in = true
    for idx in 1:3
        fn = fieldnames(TPoint)[idx]
        e1 = b.c1[idx]
        e2 = b.c2[idx]
        s = sign(e2-e1)
        is_in = is_in && (s*e1 < s*getfield(p, fn) < s*e2)
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
is_inside(x, y, z, f::Function) = z < f(x, y)
is_inside(p::TPoint, f::Function) = is_inside(p.x, p.y, p.z, f)
sample_xyz(b::Box) = rand(3).*abs.(b.c1.-b.c2).+b.minxyz

function load_spline(p)
    #=
    I don't really understand this PyNULL stuff
    I took it from here:
    https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
    Might be useful for other modules we are making
    =#
    np = PyNULL()
    copy!(np, pyimport("numpy"))
    x = np.load(p, allow_pickle=true)
    spl = x[1]
end

struct Geometry
    valley::Function
    box::Box
    tambo_center::TPoint
end

function Geometry(spl, tc::TPoint; zdown=3e4m, zup=5e3m)
    knots = spl.get_knots()
    xmin, xmax = minimum(knots[1]) * m, maximum(knots[1]) * m
    ymin, ymax = minimum(knots[2]) * m, maximum(knots[2]) * m
    # Just say that TAMBO lies on the mountain at the midpoint of the spline
    xmid, ymid = (xmax - xmin)/2, (ymax - ymin)/2
    xyzmin = [xmin-xmid, ymin-ymid, tc.z-zdown]
    xyzmax = [xmax-xmid, ymax-ymid, tc.z+zup]
    valley(x, y) = valley_helper(x, y, xyzmax.-xyzmin, spl)
    b = Box(xyzmin, xyzmax)
    Geometry(valley, b, tc)
end

function Geometry(spl; zdown=3e4m, zup=5e3m)
    knots = spl.get_knots()
    xmin, xmax = minimum(knots[1]) * m, maximum(knots[1]) * m
    ymin, ymax = minimum(knots[2]) * m, maximum(knots[2]) * m
    # Just say that TAMBO lies on the mountain at the midpoint of the spline
    xmid, ymid = (xmax - xmin)/2, (ymax - ymin)/2
    tc = TPoint(xmid, ymid, spl(xmid / m, ymid / m)[1] * m)
    xyzmin = [xmin-tc.x, ymin-tc.y, -zdown]
    xyzmax = [xmax-tc.x, ymax-tc.y, zup]
    valley(x, y) = valley_helper(x, y, tc, spl)
    b = Box(xyzmin, xyzmax)
    Geometry(valley, b, tc)
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
function valley_helper(x, y, tc, valley_spl)
    # Translate to spline coordinate system
    xt, yt = x + tc.x, y + tc.y
    # Convert to meters
    xm, ym = xt / m, yt / m
    #xm, ym = (xt |> m).val, (yt|> m).val
    # Evaluate spline and convert back to natural units
    # println(valley_spl(xm, ym)[1])
    valley_spl(xm, ym)[1]*m - tc.z
end

function (g::Geometry)(x, y)
    g.valley(x,y)[1,1]
end

function sample(n::Int, g::Geometry)
    sample(n, g.box)
end

function sample(g::Geometry)
    sample(g.box)
end

end #module