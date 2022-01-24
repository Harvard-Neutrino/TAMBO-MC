module Tracks

push!(LOAD_PATH, @__DIR__)

using Position: TPoint, Box, is_inside
export Track, normalize_track, find_position

struct Track
    ipoint::TPoint
    direction::TPoint
    function Track(ipoint::TPoint, direction::TPoint)
        new(ipoint, direction)
    end
    function Track(ipoint::TPoint, theta, phi)
        direction = TPoint(theta, phi)
        new(ipoint, direction)
    end
end

# TODO Make this geneeralize to more complex shapes
function normalize_track(t::Track, box::Box)
    if !is_inside(t.ipoint, box)
        error("Track does not start inside the box")
    end
    edges = [box.p1, box.p2]
    λ_f = Inf
    for i in 1:length(edges)
        pt = edges[i]
        prop_λs = (pt - t.ipoint)/t.direction
        for fn in fieldnames(TPoint)
            if getfield(prop_λs, fn) > 0
                λ_f = minimum((getfield(prop_λs, fn), λ_f))
            end
        end
    end
    Track(t.ipoint, t.direction*λ_f)
end

function find_position(t::Track, λ)
    t.ipoint+λ*t.direction
end

end
