struct Track{T<:Number}
    ipoint::SVector{3,T}
    fpoint::SVector{3,T}
    direction::Direction
    norm::Float64
    function Track(ipoint::T, fpoint::U, d::Direction, norm::Float64) where {T,U}
        ipoint, fpoint = promote(ipoint, fpoint)
        return new{eltype(ipoint)}(ipoint, fpoint, d, norm)
    end
end

function Track(fpoint::SVector{3})
    ipoint = SVector{3}([0, 0, 0])
    d = Direction(fpoint, ipoint)
    l = norm(fpoint .- ipoint)
    return Track(ipoint, fpoint, d, l)
end

function Track(ipoint::SVector{3,T}, fpoint::SVector{3,U}) where {T,U}
    promote(ipoint, fpoint)
    d = Direction(fpoint, ipoint)
    l = norm(fpoint .- ipoint)
    return Track(ipoint, fpoint, d, l)
end

function Track(ipoint::SVector{3}, d::Direction, b::Box)
    fpoint = intersect(ipoint, d, b)
    l = norm(fpoint .- ipoint)
    return Track(ipoint, fpoint, d, l)
end

function (t::Track)(λ)
    return t.ipoint .+ λ * t.norm * t.direction
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
    return λf * d .+ p
end

function Base.intersect(t::Track, b::Box)
    return intersect(t.ipoint, t.direction, b)
end

function Base.intersect(t::Track, z::Float64)
    Δz = t.fpoint.z - t.ipoint.z
    root = (z - t.ipoint.z) / Δz
    if root < 0
        root = Inf
    elseif root > 1
        root = Inf
    end
    return root
end

function Base.intersect(t::Track, v::Function)
    oned_valley = reduce_f(t, v)
    root_func(λ) = oned_valley(λ) - t(λ).z
    return zeros = find_zeros(root_func, 0, 1)
end

function Base.intersect(t::Track, g::Geometry)
    ixs = intersect.(Ref(t), g.zboundaries)
    ixs = vcat(ixs, intersect(t, g.valley))
    ixs = ixs[ixs .<= 1]
    return ixs = sort(ixs)
end

function Base.reverse(t::Track)
    return Track(t.fpoint, t.ipoint)
end

function reduce_f(t::Track, f)
    return g(λ) = f(t(λ).x, t(λ).y)
end

####### Range struct #######

struct Range
    λstart::Float64
    λwidth::Float64
    pstart::SVector{3}
    length::Float64
    density::Float64
    medium_name::String
end

function Range(λstart::Float64, λwidth::Float64, t::Track, density::Float64)
    pstart = t(λstart)
    width = t.norm * λwidth
    # I'm sorry whoever finds this....
    medium_name = density > 1 ? "StandardRock" : "Air"
    range = Range(λstart, λwidth, pstart, width, density, medium_name)
    return range
end

function Base.show(io::IO, r::Range)
    return print(
        io,
        """
        λstart: $(r.λstart)
        λwidth: $(r.λwidth)
        pstart (m): $(r.pstart / units.m)
        width (m): $(r.length / units.m)
        density (g/cm^3): $(r.density / (units.gr / units.cm^3))
        medium_name: $(r.medium_name)
        """,
    )
end

function inversecolumndepth(tr::Track, cd::T, g::Geometry) where {T<:Number}
    ranges = computeranges(tr, g)
    return inversecolumndepth(tr, cs, valley, ranges)
end

function inversecolumndepth(tr::Track, cd::T, g::Geometry, ranges::Vector) where {T<:Number}
    f(λ) = columndepth(tr, λ, ranges) - cd
    λ_int = find_zero(f, (0, 1))
    return λ_int
end

function computeranges(t::Track, g::Geometry, ixs)
    rgen = vcat(0, ixs, 1)
    ranges = [
        Range(x[1], x[2] - x[1], t, getdensity(t((x[1] + x[2]) / 2), g)) for
        x in zip(rgen[1:(end - 1)], rgen[2:end])
    ]
    return ranges
end

function computeranges(t::Track, g::Geometry)
    ixs = intersect(t, g)
    return computeranges(t, g, ixs)
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
    return ρ
end

function columndepth(t::Track, λ::T, ranges::Vector) where {T<:Number}
    ranges = [r for r in ranges if r.λstart < λ]
    if length(ranges) > 0
        endrange = ranges[end]
        ranges[end] = Range(endrange.λstart, λ - endrange.λstart, t, endrange.density)
    end
    cd = 0units[:mwe]
    for r in ranges
        cd += r.length * r.density
    end
    return cd
end

function columndepth(t::Track, λ::T, g::Geometry) where {T<:Number}
    ranges = computeranges(t, g)
    return columndepth(t, λ, ranges)
end

function totalcolumndepth(t::Track, ranges::Vector)
    return columndepth(t, 1, ranges)
end

function totalcolumndepth(t::Track, g::Geometry)
    ranges = computeranges(t, g)
    return totalcolumndepth(t, ranges)
end
