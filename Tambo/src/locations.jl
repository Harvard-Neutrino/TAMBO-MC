const testsite_coord = Coord(deg2rad(-15.58714), deg2rad(-71.9765237))
const minesite_coord = Coord(deg2rad(-15.664653), deg2rad(-72.1547479))
const minesite_normal_vec = Direction([0.732001, 0.59897, 0.324666])

struct Coord{T<:Float64}
    latitude::T
    longitude::T
end

function Base.show(io::IO, coord::Coord)
    roundlat = round(coord.latitude*180/π, sigdigits=5)
    roundlong = round(coord.longitude*180/π, sigdigits=5)
    print(
        io,
        """
        ($(roundlat)°, $(roundlong)°)"""
    )
end

function mincoord_fromfile(filename)
    open(filename) do file
        latmin = 0
        longmin = 0
        for (idx, line) in enumerate(eachline(file))
            if idx!=1
                splitline = split(line, "\t")
                lat = parse(Float64, splitline[2])
                long = parse(Float64, splitline[3])
                latmin = minimum([lat, latmin])
                longmin = minimum([long, longmin])
            end
        end
        return Coord(deg2rad(latmin), deg2rad(longmin))
    end
end

function mincoord_fromfile()
    filename = "$(@__DIR__)/../resources/ColcaValleyData.txt"
    return mincoord_fromfile(filename)
end

"""
    latlong_to_xy(lat, long, latmin, longmin)

TBW
"""
function latlong_to_xy(lat, long, latmin, longmin)
    r = 6_371.0 * units.km
    x = r * cos(latmin) * (long - longmin)
    y = r * (lat - latmin)
    return SVector{2}(x, y)
end

"""
    latlong_to_xy(coord::Coord, coordmin::Coord)

TBW
"""
function latlong_to_xy(coord::Coord, coordmin::Coord)
    return latlong_to_xy(
        coord.latitude,
        coord.longitude,
        coordmin.latitude,
        coordmin.longitude
    )
end

"""
    xy_to_latlong(x, y, latmin, longmin)

TBW
"""
function xy_to_latlong(x, y, latmin, longmin)
    r = 6_371.0 * units.km
    long = x / (r * cos(latmin)) + longmin
    lat = y / r + latmin
    return Coord(lat, long)
end

"""
    xy_to_latlong(xy, latmin, longmin)

TBW
"""
function xy_to_latlong(xy, latmin, longmin)
    return xy_to_latlong(xy[1], xy[2], latmin, longmin)
end

"""
    xy_to_latlong(xy, coordmin::Coord)

TBW
"""
function xy_to_latlong(xy, coordmin::Coord)
    return xy_to_latlong(xy, coordmin.latitude, coordmin.longitude)
end