struct Coord{T<:Float64}
    latitude::T
    longitude::T
end

function mincoord_fromfile(;filename = "$(@__DIR__)/../resources/ColcaValleyData.txt")
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

const testsite_coord = Coord(deg2rad(-15.58714), deg2rad(-71.9765237))
const minesite_coord = Coord(deg2rad(-15.664653), deg2rad(-72.1547479))

"""
    latlong_to_xy(lat, long, latmin, longmin)

Function to calculate the xy coordinate for a latitude and longitude
coordinate. This implicitly assumes that Δθ and Δϕ are small so that
the sphere is locally floating. All angles are in radians !!!!!
latmin and longmin are -15.73975004° and -72.336236836° for the current spline

# Example
```julia-repl
julia> latmin, longmin = deg2rad(-15.73975004), deg2rad(-72.336236836);

julia> tapay_lat, tapay_long = deg2rad.((-15.63, -72.16));

julia> x_tapay, y_tapay = latlong_to_xy(tapay_lat, tapay_long, latmin, longmin) ./ units.km
(18.861842060081564, 12.203647647037522)
```
"""
function latlong_to_xy(lat, long, latmin, longmin)
    r = 6_371.0 * units.km
    x = r * cos(latmin) * (long - longmin)
    y = r * (lat - latmin)
    return SVector{2}(x, y)
end
#function latlong_to_xy(lat, long, latmin, longmin)
#    r = 6_371.0 * units.km
#    x = r * cos(latmin) * (long - longmin)
#    y = r * (lat - latmin)
#    return SVector{2}([x, y])
#end

function latlong_to_xy(coord::Coord, coordmin::Coord)
    return latlong_to_xy(
        coord.latitude,
        coord.longitude,
        coordmin.latitude,
        coordmin.longitude
    )
end