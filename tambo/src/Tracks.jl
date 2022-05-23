module Tracks

push!(LOAD_PATH, @__DIR__)

using Geometries
using Units
using Unitful
using LinearAlgebra
using Roots: find_zeros, find_zero
export Track, total_column_depth, column_depth, inverse_column_depth, intersect


struct Track
    ipoint::TPoint
    fpoint::TPoint
    direction::Direction
    norm
end

function Track(fpoint::TPoint)
    ipoint = TPoint(0, 0, 0)
    d = Direction(fpoint, ipoint)
    l = norm(fpoint-ipoint)
    Track(ipoint, fpoint, d, l)
end

function Track(ipoint::TPoint, fpoint::TPoint)
    d = Direction(fpoint, ipoint)
    l = norm(fpoint-ipoint)
    Track(ipoint, fpoint, d, l)
end

function Track(ipoint::TPoint, d::Direction, b::Box)
    fpoint = intersect(ipoint, d, b)
    l = norm(fpoint-ipoint)
    Track(ipoint, fpoint, d, l)
end

function (t::Track)(λ)
    t.ipoint+λ*t.norm*t.direction
end

function Base.intersect(t::Track, b::Box)
    intersect(t.ipoint, t.direction, b)
end
∩(t::Track, b::Box) = intersect(t, b)

function Base.intersect(p::TPoint, d::Direction, box::Box)
    if !is_inside(p, box)
        error("Track does not start inside the box")
    end
    edges = [box.c1, box.c2]
    #=
    λ is proportional to the distance the particle propagates
    We will use distance interchangably with λ
    =#
    λ_f = Inf*m
    # Iterate over points which define box
    for edge in edges
        # Points where track shares coordinate with box edges
        prop_λs = (edge - p)/d
        # Iterate over x, y, z
        for fn in fieldnames(TPoint)
            λ = getfield(prop_λs, fn)
            #= 
            λ > 0 means the track is moving forward
            λ < 0 means the track is moving backwards
            λ = 0 means the tracks starts on an edge of the box
            We want a forward-going track
            =#
            if λ >= 0m
                #=
                The track hits the box at minimum prop distance
                Else it will be outside the box before travelling that distance
                =#
                λ_f = minimum((λ, λ_f))
            end
        end
    end
    λ_f*d+p
end
∩(p::TPoint, d::Direction, b::Box) = intersect(p, d, b)

function Base.reverse(t::Track)
    Track(t.fpoint, t.ipoint)
end

function reduce_f(t::Track, f)
    g(λ) = f(t(λ).x, t(λ).y)
end

function inverse_column_depth(
    tr::Track, cd, valley; ρ_air=ρ_air0, ρ_rock=ρ_rock0, ixs=Nothing
)
    f(λ) = column_depth(tr, λ, valley; ρ_air=ρ_air0, ρ_rock=ρ_rock0, ixs=ixs) - cd
    λ_int = find_zero(f, (0,1))
    λ_int
end

function compute_ranges(t, valley, ixs, λ)
    ixs = ixs[ixs .< λ]
    rgen = vcat([0], ixs, [λ])
    ranges = [
        (x[2]-x[1], is_inside(t((x[1]+x[2])/2), valley))
        for x in zip(rgen[1:end-1], rgen[2:end])
    ]
end

function Base.intersect(t::Track, v)
    oned_valley = reduce_f(t, v)
    root_func(λ) = oned_valley(λ)-t(λ).z
    zeros = find_zeros(root_func, 0, 1)
end
#∩(t::Track, v) = intersect(t, v)

function column_depth(t, λ, valley; ρ_air=ρ_air0, ρ_rock=ρ_rock0, ixs=Nothing)
    if ixs==Nothing
        ixs = intersect(t, valley)
    end
    ranges = compute_ranges(t, valley, ixs, λ)
    cd = 0mwe
    for x in ranges
        width, in_mountain = x
        ρ = in_mountain ? ρ_rock : ρ_air
        cd += width*t.norm*ρ
    end
    cd
end

function total_column_depth(t, valley; ρ_air=ρ_air0, ρ_rock=ρ_rock0, ixs=Nothing)
    column_depth(t, 1, valley, ixs=ixs)
end

end # module
