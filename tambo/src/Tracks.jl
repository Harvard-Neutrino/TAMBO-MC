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
    function Track(ipoint::TPoint, theta::AbstractFloat, phi::AbstractFloat)
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
    #=
    λ is proportional to the distance the particle propagates
    We will use distance interchangably with λ
    =#
    λ_f = Inf
    # Iterate over points which define box
    for i in 1:length(edges)
        pt = edges[i]
        # Points where track shares coordinate with box edges
        prop_λs = (pt - t.ipoint)/t.direction
        # Iterate over x, y, z
        for fn in fieldnames(TPoint)
            #= 
            λ > 0 means the track is moving forward
            λ < 0 means the track is moving backwards
            λ = 0 means the tracks starts on an edge of the box
            We want a forward-going track
            =#
            if getfield(prop_λs, fn) >= 0
                #=
                The track hits the box at minimum prop distance
                Else it will be outside the box before travelling that distance
                =#
                λ_f = minimum((getfield(prop_λs, fn), λ_f))
            end
        end
    end
    Track(t.ipoint, t.direction*λ_f)
end

function complement(t::Track, b::Box)
    nt = normalize_track(t, b)
    # Make a track which starts on the edge of the box and points in
    Track(find_position(nt, 1), -1*t.direction)
end

function find_position(t::Track, λ)
    t.ipoint+λ*t.direction
end

end
