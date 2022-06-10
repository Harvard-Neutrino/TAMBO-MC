include("geometries.jl")
include("units.jl")

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

function Base.intersect(p::SVector{3}, d::Direction, box::Box)
    if !inside(p, box)
        error("Track does not start inside the box")
    end
    edges = [box.c1, box.c2]
    λf = Inf
    # Iterate over points which define box
    for edge in edges
        # Points where track shares coordinate with box edges
        prop_λs = Vector((edge .- p) ./ d.proj)
        # Negative λ are backwards moving which is not what we want
        prop_λs[prop_λs .<= 0] .= Inf
        λf = minimum((λf, minimum(prop_λs)))
    end
    λf * d .+ p
end

function Base.intersect(t::Track, b::Box)
    intersect(t.ipoint, t.direction, b)
end

function Base.intersect(t::Tambo.Track, z::Float64)
    Δz = t.fpoint.z - t.ipoint.z
    root = (z-t.ipoint.z) / Δz
    if root < 0
        root = Inf
    elseif root > 1
        root = Inf
    end
    root
end

function Base.intersect(t::Track, v::Function)
    oned_valley = reduce_f(t, v)
    root_func(λ) = oned_valley(λ)-t(λ).z
    zeros = find_zeros(root_func, 0, 1)
end

function Base.intersect(t::Track, g::Geometry)
    ixs = intersect.(Ref(t), g.zboundaries)
    ixs = vcat(ixs, intersect(t, g.valley))
    ixs = ixs[ixs.<=1]
    ixs = sort(ixs)
end

function Base.reverse(t::Track)
    Track(t.fpoint, t.ipoint)
end

function reduce_f(t::Track, f)
    g(λ) = f(t(λ).x, t(λ).y)
end

function inversecolumndepth(tr::Track, cd::T, g::Geometry) where {T<:Number}
    ranges = computeranges(tr, g)
    inversecolumndepth(tr, cs, valley, ranges)
end

function inversecolumndepth(tr::Track, cd::T, g::Geometry, ranges::Vector) where {T<:Number}
    f(λ) = columndepth(tr, λ, ranges) - cd
    λ_int = find_zero(f, (0,1))
    λ_int
end

function computeranges(t::Track, g::Geometry, ixs)
    rgen = vcat(0, ixs, 1)
    ranges = [
        (x[1], x[2]-x[1], getdensity(t((x[1]+x[2])/2), g))
        for x in zip(rgen[1:end-1], rgen[2:end])
    ]
end

function computeranges(t::Track, g::Geometry)
    ixs = intersect(t, g)
    computeranges(t, g, ixs)
end

function getdensity(p::SVector{3}, g::Geometry)
    ρ = 0
    if !inside(p, g.valley)
        ρ = g.ρair
    elseif any(p.z .< g.zboundaries)
        ρ = g.ρs[findlast(p.z .< g.zboundaries)]
    else
        ρ = g.ρrock
    end
    ρ
end

function columndepth(t::Track, λ::T, ranges::Vector) where {T<:Number}
    ranges = [r for r in ranges if r[1]<λ]
    if length(ranges)>0
        ranges[end] = (ranges[end][1], λ - ranges[end][1], ranges[end][3])
    end
    cd = 0units[:mwe]
    for (_, width, ρ) in ranges
        cd += width * t.norm * ρ
    end
    cd
end

function columndepth(t::Track, λ::T, g::Geometry) where {T<:Number}
    ranges = computeranges(t, g)
    columndepth(t, λ, ranges)
end

function totalcolumndepth(t::Track, ranges::Vector)
    columndepth(t, 1, ranges)
end

function totalcolumndepth(t::Track, g::Geometry)
    ranges = computeranges(t, g)
    totalcolumndepth(t, ranges)
end