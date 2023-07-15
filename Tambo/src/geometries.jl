const defaults = (zup=5units.km, zdown=50units.km, depth1=12units.km, depth2=21.4units.km)

struct Plane
    n̂::Direction
    x0::SVector{3, Float64}
    function Plane(n, x0)
        n̂ = Direction(n)
        x0 = SVector{3, Float64}(x0)
        return new(n̂, x0)
    end
end

function Base.show(io::IO, plane::Plane)
    s = "Plane(n̂=$(plane.n̂), x0=$(plane.x0 / units.m))"
    print(io, s)
end

struct Box
    c1::SVector{3, Float64}
    c2::SVector{3, Float64}
    function Box(c1, c2)
        c1, c2 = SVector{3, Float64}(c1), SVector{3, Float64}(c2)
        return new(c1, c2)
    end
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

struct Valley
    spline::Spline2D
    mincoord::Coord
end

function Valley(splpath::String)
    f = jldopen(splpath)
    spl = f["spline"]
    mincoord = f["mincoord"]
    return Valley(spl, mincoord)
end    



#function load_spline(p, key="spline")
#    f = jldopen(p, "r")
#    spl = f[key]
#    return spl
#end

struct Geometry
    valley::Valley
    box::Box
    tambo_offset::SVector{3, Float64}
    ρair::Float64
    ρrock::Float64
    tambo_bounds::SMatrix{4, 2, Float64}
end

function Geometry(spl_path::String, tambo_coord::Coord)
    spl = nothing
    min_coord = nothing
    jldopen(spl_path) do jldf
        spl = jldf["spline"]
        min_coord = jldf["mincoord"]
    end
    valley = Valley(spl, min_coord)
    tambo_xy = latlong_to_xy(tambo_coord, min_coord)

    z = spl(tambo_xy.x / units.m, tambo_xy.y / units.m) * units.m
    xyzoffset = SVector{3}([tambo_xy.x, tambo_xy.y, z])
    zup = defaults.zup
    zdown = defaults.zdown
    knots = spl.tx, spl.ty
    xmin, xmax = minimum(knots[1]) * units.m, maximum(knots[1]) * units.m
    ymin, ymax = minimum(knots[2]) * units.m, maximum(knots[2]) * units.m
    xyzmin = SVector{3}([xmin - xyzoffset.x, ymin - xyzoffset.y, -zdown])
    xyzmax = SVector{3}([xmax - xyzoffset.x, ymax - xyzoffset.y, zup])
    box = Box(xyzmin, xyzmax)
    bounds = units.km .* SMatrix{4,2, Float64}([
        -1.0 1.0 -1.0 1.0;
        -0.5 -0.5 0.5 0.5
    ])
    return Geometry(valley, box, xyzoffset, units.ρair0, units.ρrock0, bounds)
end

function Plane(n̂::Direction, coord::Coord, geo::Geometry)
    plane_xy = latlong_to_xy(coord, geo.valley.mincoord) .- geo.tambo_offset[1:2]
    plane_z = geo(plane_xy...)
    x0 = SVector{3}(
        plane_xy...,
        plane_z
    )
    return Plane(n̂.proj, x0)
end

function (geo::Geometry)(x, y)
    return valley_helper(x, y, geo.tambo_offset, geo.valley.spline)
end

function (geo::Geometry)(sv::SVector{2})
    return geo(sv...)
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
```
"""
function valley_helper(x, y, xyzoffset, valley_spl)
    # Translate to spline coordinate system
    xt, yt = x + xyzoffset[1], y + xyzoffset[2]
    # Convert to meters
    xm, ym = xt / units.m, yt / units.m
    # Evaluate spline and convert back to natural units
    return valley_spl(xm, ym) * units.m - xyzoffset[3]
end

"""
    density(p::SVector{3}, g::Geometry)

TBW
"""
function density(p::SVector{3}, geo::Geometry)
    ρ = 0
    if !inside(p, geo)
        ρ = geo.ρair
    else
        ρ = geo.ρrock
    end
    return ρ
end

function inside(sv::SVector{3}, b::Box)
    e1 = b.c1
    e2 = b.c2
    s = sign.(e2 .- e1)
    is_in = all(s .* e1 .< s .* sv < s .* e2)
    return is_in
end
inside(x, y, z, f::Function) = z < f(x, y)
inside(x, y, z, geo::Geometry) = z < geo(x, y)
inside(sv::SVector{3}, f::Function) = inside(sv.x, sv.y, sv[3], f)
inside(sv::SVector{3}, geo::Geometry) = inside(sv.x, sv.y, sv[3], geo)
inside(p::Particle, f::Function) = inside(p.position, f)
inside(ps::Vector{Particle}, f::Function) = [inside(p.position, f) for p in ps]