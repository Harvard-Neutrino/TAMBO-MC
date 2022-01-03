
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

struct Point{T<:Real}
    longitude::T
    latitude::T
    x::T
    y::T
    z::T
    function Point(longitude::T, latitude::T, elevation::T; latmin = 0.0, longmin = 0.0) where T<:Real
    x,y = get_distance_coordinate(latitude,longitude,latmin, longmin)
    new{T}(longitude, latitude, x, y, elevation)
    end
end

struct Direction{T<:Real}
    phi::T
    theta::T
    x::T
    y::T
    z::T
    function Direction{T<:Real}(phi::T, theta::T)
        x = cos(theta)*sin(phi)
        y = sin(theta)*sin(phi)
        z = cos(phi)
        new{T}(phi,theta,x,y,z)
    end
end

# original FORTRAN interpolation library
# that is inside scipy interpolate
using Dierckx
#using CSV,DataFrames,LinearAlgebra
using DelimitedFiles

struct Geometry
    lat_min::Float
    lat_max::Float
    long_min::Float
    long_max::Float
    elev_min::Float
    elev_max::Float
    geometry_spline::Spline2D
    geometry_box
    density_rock::Float #"kg/m^3"
    density_air::Float #"kg/m^3" 
    function Geometry(textfile)
        data_input, header_input = readdlm("./ColcaValleyData.txt",'\t', '\n', header = true)
        latitude = deg2rad(data_input[:,2])
        longitude = deg2rad(data_input[:,3])
        elevation = data_input[:,4]

        lat_min = minimum(latitude)
        lat_max = maximum(latitude)
        long_min = minimum(longitude)
        long_max = maximum(longitude)
        elev_min = minimum(elevation)
        elev_max = maximum(elevation)

        number_geo_points = len(latitude)

        # A coord. system in meters with the origin at latmin, latmax
        coordinate_points = ([Point(longitude[i],
                                    latitude[i],
                                    elevation[i];
                                    lat_min, 
                                    long_min) 
        for i in 1:number_geo_points])

        geometry_spline = Spline2D(map(p->p.x,coordinate_points),
                                   map(p->p.y,coordinate_points),
                                   map(p->p.z,coordinate_points);
                                   kx=3,ky=3
                                   )

        x_max = max(map(p->p.x,coordinate_points))
        y_max = max(map(p->p.y,coordinate_points))
        geometry_box = [0.0,x_max,0.0,y_max,elev_min,elev_max]

        density_rock =  5520 #"kg/m^3"
        density_air = 1225 #"kg/m^3" 

        new(lat_min,lat_max,long_min,long_max,elev_min,elev_max,geometry_spline,geometry_box,density_rock,density_air)
    end
end