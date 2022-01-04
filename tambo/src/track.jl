include("./geometry.jl")
using .geometry


function get_start_track(point, dir)
    #=
    Builds 6 item array defining starting position and direction of track
    =#    
    line_eq = [[point.x, point.y, point.z],
               [dir.x, dir.y, dir.z]]
    round.(line_eq, 3)
end

function find_end_points(t::Track, geo::Geometry)
    
end

struct Track
    ipoint::Point
    direction::Direction
end

