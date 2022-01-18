module track

    include("./position.jl")
    using .position
    export Track, TPoint, is_inside, normalize_track, find_position

    function is_inside(p::TPoint, box::Vector{TPoint})
        is_it = true
        for fn in fieldnames(TPoint)
            edge_1 = getfield(box[1], fn)
            edge_2 = getfield(box[2], fn)
            s = sign((edge_2-edge_1))
            is_it = is_it && (s*edge_1 < s*getfield(p, fn) < s*edge_2)
            if !is_it
                return false
            end
        end
        is_it
    end

    struct Track
        ipoint::TPoint
        direction::TPoint
    end

    function normalize_track(t::Track, box::Vector{TPoint})
        if !is_inside(t.ipoint, box)
            error("Track does not start inside the box")
        end
        λ_f = Inf
        for i in 1:length(box)
            pt      = box[i]
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
