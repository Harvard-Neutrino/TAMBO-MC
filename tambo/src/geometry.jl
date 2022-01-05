module geometry

export Geometry

# ----------------------------------------------------------------------------
# original FORTRAN interpolation library
# that is inside scipy interpolate
using Dierckx
#using CSV,DataFrames,LinearAlgebra
using DelimitedFiles

struct Geometry
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
    function Geometry(textfile::String)
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

end #module

if abspath(PROGRAM_FILE) == @__FILE__
    geo = geometry.Geometry("../../resources/ColcaValleyData.txt")
    println("hola")
end