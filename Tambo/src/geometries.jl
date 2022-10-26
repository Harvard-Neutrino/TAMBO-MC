const defaults = (zup=5units.km, zdown=50units.km, depth1=12units.km, depth2=21.4units.km)

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