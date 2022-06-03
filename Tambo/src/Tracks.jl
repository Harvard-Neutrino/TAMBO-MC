include("Geometries.jl")
include("Units.jl")

using LinearAlgebra: norm
using Roots: find_zeros, find_zero
using StaticArrays

struct Track{T<:Number}
    ipoint::SVector{3, T}
    fpoint::SVector{3, T}
    direction::Direction
    norm::Float64
    function Track(ipoint::T, fpoint::U, d::Direction, norm::Float64) where {T,U}
        ipoint, fpoint = promote(ipoint, fpoint)
        new{eltype(ipoint)}(ipoint, fpoint, d, norm)
    end
end

function Track(fpoint::SVector{3})
    ipoint = SVector{3}([0, 0, 0])
    d = Direction(fpoint, ipoint)
    l = norm(fpoint .- ipoint)
    Track(ipoint, fpoint, d, l)
end

function Track(ipoint::SVector{3, T}, fpoint::SVector{3, U}) where {T, U}
    promote(ipoint, fpoint)
    d = Direction(fpoint, ipoint)
    l = norm(fpoint .- ipoint)
    Track(ipoint, fpoint, d, l)
end

function Track(ipoint::SVector{3}, d::Direction, b::Box)
    fpoint = intersect(ipoint, d, b)
    l = norm(fpoint .- ipoint)
    Track(ipoint, fpoint, d, l)
end

function (t::Track)(λ)
    t.ipoint .+ λ*t.norm*t.direction
end

function Base.intersect(t::Track, b::Box)
    intersect(t.ipoint, t.direction, b)
end

function Base.intersect(p::SVector{3}, d::Direction, box::Box)
    if !is_inside(p, box)
        error("Track does not start inside the box")
    end
    edges = [box.c1, box.c2]
    #=
    λ is proportional to the distance the particle propagates
    We will use distance interchangably with λ
    =#
    λf = Inf*units[:m]
    # Iterate over points which define box
    for edge in edges
        # Points where track shares coordinate with box edges
        # TODO Make ditection have a SVector and fix this
        prop_λs = Vector((edge .- p) ./ d.proj)
        prop_λs[prop_λs .<= 0] .= Inf
        #= 
        λ > 0 means the track is moving forward
        λ < 0 means the track is moving backwards
        λ = 0 means the tracks starts on an edge of the box
        We want a forward-going track
        =#
        λf = minimum((λf, minimum(prop_λs)))
    end
    λf * d .+ p
end
∩(p::SVector{3}, d::Direction, b::Box) = intersect(p, d, b)

function Base.reverse(t::Track)
    Track(t.fpoint, t.ipoint)
end

function reduce_f(t::Track, f)
    g(λ) = f(t(λ).x, t(λ).y)
end

function inverse_column_depth(
    tr::Track,
    cd,
    valley;
    ρ_air=units[:ρ_air0],
    ρ_rock=units[:ρ_rock0]
)
    ixs = intersect(tr, valley)
    inverse_column_depth(tr, cs, valley, ixs)
end

function inverse_column_depth(
    tr::Track,
    cd,
    valley,
    ixs;
    ρ_air=units[:ρ_air0],
    ρ_rock=units[:ρ_rock0]
)
    f(λ) = column_depth(tr, λ, valley, ixs; ρ_air=ρ_air, ρ_rock=ρ_rock) - cd
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


function column_depth(
    t,
    λ,
    valley,
    ixs;
    ρ_air=units[:ρ_air0],
    ρ_rock=units[:ρ_rock0],
)
    ranges = compute_ranges(t, valley, ixs, λ)
    cd = 0units[:mwe]
    for x in ranges
        width, in_mountain = x
        ρ = in_mountain ? ρ_rock : ρ_air
        cd += width*t.norm*ρ
    end
    cd
end

function column_depth(
    t,
    λ,
    valley;
    ρ_air=units[:ρ_air0],
    ρ_rock=units[:ρ_rock0],
)
    ixs = intersect(t, valley)
    column_depth(t, λ, valley, ixs)
end

function total_column_depth(
    t,
    valley,
    ixs;
    ρ_air=units[:ρ_air0],
    ρ_rock=units[:ρ_rock0],
)
    column_depth(t, 1.0, valley, ixs; ρ_air=ρ_air, ρ_rock=ρ_rock)
end

function total_column_depth(
    t,
    valley;
    ρ_air=units[:ρ_air0],
    ρ_rock=units[:ρ_rock0],
)
    ixs = intersect(t, valley)
    total_column_depth(t, valley, ixs; ρ_air=ρ_air, ρ_rock=ρ_rock)
end