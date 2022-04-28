module Geometry

using LinearAlgebra, StaticArrays
using PyCall
using Units: m

export TAMBOGeometry,
       Coord, 
       TPoint,
       Box,
       is_inside,
       GenerationRegion,
       sample_xyz

# ----------------------------------------------------------------------------
# original FORTRAN interpolation library
# that is inside scipy interpolate
using Dierckx
#using CSV,DataFrames,LinearAlgebra
using DelimitedFiles

struct TAMBOGeometry
    lat_min::Float64
    lat_max::Float64
    long_min::Float64
    long_max::Float64
    elev_min::Float64
    elev_max::Float64
    geometry_spline::Spline2D
    geometry_box::Vector{Float64}
    density_rock::Float64 #"kg/m^3"
    density_air::Float64 #"kg/m^3" 
    function TAMBOGeometry(textfile::String)
        data_input, header_input = readdlm(textfile,'\t', '\n', header = true)
        latitude = map(deg2rad,data_input[:,2])
        longitude = map(deg2rad,data_input[:,3])
        elevation = data_input[:,4]

        lat_min = minimum(latitude)
        lat_max = maximum(latitude)
        long_min = minimum(longitude)
        long_max = maximum(longitude)
        elev_min = minimum(elevation)
        elev_max = maximum(elevation)

        number_geo_points = length(latitude)

        # A coord. system in meters with the origin at latmin, latmax
        coordinate_points = [TPoint(longitude[i],
                                    latitude[i],
                                    elevation[i];
                                    lat_min, 
                                    long_min) 
        for i in 1:number_geo_points]
        geometry_spline = Spline2D(map(p->p.x,coordinate_points),
                                   map(p->p.y,coordinate_points),
                                   map(p->p.z,coordinate_points);
                                   kx=3,ky=3,s=100_000
                                   )

        x_max = max(map(p->p.x,coordinate_points))
        y_max = max(map(p->p.y,coordinate_points))
        geometry_box = [0.0,x_max,0.0,y_max,elev_min,elev_max]

        density_rock = 5520 #"kg/m^3"
        density_air = 1225 #"kg/m^3" 

        new(lat_min,lat_max,long_min,long_max,elev_min,elev_max,geometry_spline,geometry_box,density_rock,density_air)
    end
end



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
#Base.:*(t, p::TPoint) = p*t
LinearAlgebra.norm(p::TPoint) = norm((p.x, p.y, p.z))
#Base.:*(r::Rotations.Rotation, p::TPoint) = TPoint(r*[p.x, p.y, p.z]...)
#function (t::Translation)(p::TPoint)
#    TPoint(t([p.x, p.y, p.z])...)
#end


struct Box
    c1::SVector{3}
    c2::SVector{3}
    minxyz::SVector{3}
    maxxyz::SVector{3}
    function Box(c1, c2)
        c1 = SVector{3}(c1)
        c2 = SVector{3}(c2)
        minxyz = SVector{3}(minimum.(zip(c1,c2)))
        maxxyz = SVector{3}(maximum.(zip(c1,c2)))
        new(c1, c2, minxyz, maxxyz)
    end
    function Box(x, y, z)
        c1 = SVector{3}([0,0,0])
        c2 = SVector{3}([x,y,z])
        new(c1, c2)
    end
    function Box(c)
        c1 = SVector{3}([0,0,0])
        c2 = SVector{3}(c)
        new(c1, c2)
    end
end

function sample(n::Int, b::Box)
    [rand(3).*abs.(b.c1.-b.c2).+b.minxyz for _ in 1:n]
end

function sample(b::Box)
    rand(3).*abs.(b.c1.-b.c2).+b.minxyz
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

#=
This does not work yet because we do not have a load
=#
struct GenerationRegion
    valley::Function
    box::Box
    function GenerationRegion(spline)
        new(construct_generation_region(spline)...)
    end
    function GenerationRegion(spline_path::String)
        # Load the spline first
        spline = load_spline(spline_path)
        new(construct_generation_region(spline)...)
    end
end

function (gr::GenerationRegion)(x, y)
    gr.valley(x,y)[1,1]
end

"""
    valley_helper(x, y)

Function for calling the Python `valley_spl` using unitful Julia 
`x` and `y` coordinates. It also accounts for the coordinate change `Δc` between
the spline coordinates and TAMBO coordinates

# Example
```julia-repl
julia> valley_helper(1km, 4000m, [1000m, 1000m, 1km], spl)
```
"""
function valley_helper(x, y, Δc, valley_spl)
    xt, yt = x + Δc[1]/2, y + Δc[2]
    xm, ym = xt / m, yt / m
    #xm, ym = (xt |> m).val, (yt|> m).val
    valley_spl(xm, ym)[1]*m + Δc[3]
end

function construct_generation_region(spl, zmin=-3e4m, zmax=5e3m)
    #=
    This function will need to return a spline which takes an x and y coordinate 
    and returns the height of the valley at that point.
    It will also need to return a Box which defines the generation region
    valley, box
    =#
    knots = spl.get_knots()
    xmin, xmax = minimum(knots[1]), maximum(knots[1])
    ymin, ymax = minimum(knots[2]), maximum(knots[2])
    #=
    I think that this only works if the axes of our box
    are aligned in a certain way. We may want to do something more general
    =#
    cmin = [xmin*m, ymin*m, zmin]
    cmax = [xmax*m, ymax*m, zmax]
    cmid = (cmin .+ cmax)/2
    Δc = cmax .- cmin
    # Make a box which is the size of the spline region and is centered at 0
    b = Box(-Δc/2, Δc/2)
    valley(x, y) = valley_helper(x, y, Δc, spl)
    valley, b
end


function sample(n::Int, gr::GenerationRegion)
    sample(n, gr.box)
end

function sample(gr::GenerationRegion)
    sample(gr.box)
end

end #module

if abspath(PROGRAM_FILE) == @__FILE__
    geo = geometry.TAMBOGeometry("../../resources/ColcaValleyData.txt")
end