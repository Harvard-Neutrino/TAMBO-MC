module Tracks

push!(LOAD_PATH, @__DIR__)

using Geometry: TPoint, Box, is_inside
using Units
using Unitful
using LinearAlgebra
using Roots: find_zeros, find_zero
export Track, Direction, total_column_depth, column_depth, inverse_column_depth

struct Direction
    θ
    ϕ
    x_proj
    y_proj
    z_proj
    function Direction(θ, ϕ)
        θ, ϕ = Float64(θ), Float64(ϕ)
        x = cos(ϕ)*sin(θ)
        y = sin(ϕ)*sin(θ)
        z = cos(θ)
        new(θ, ϕ, x, y, z)
    end
    function Direction(x, y, z)
        x,y,z = (x,y,z)./norm((x,y,z))
        θ = acos(z)
        ϕ = atan(y, x)
        new(θ, ϕ, x, y, z)
    end
    function Direction(p::TPoint)
        Direction(p.x, p.y, p.z)
    end
    function Direction(ip::TPoint, fp::TPoint)
        Direction(ip-fp)
    end
end

Base.:/(p::TPoint, d::Direction) = TPoint(p.x/d.x_proj, p.y/d.y_proj, p.z/d.z_proj)
Base.:*(m, d::Direction) = TPoint(m*d.x_proj, m*d.y_proj, m*d.z_proj)

struct Track
    ipoint::TPoint
    fpoint::TPoint
    direction::Direction
    norm
    function Track(ipoint::TPoint, fpoint::TPoint)
        d = Direction(fpoint, ipoint)
        l = norm(fpoint-ipoint)
        new(ipoint, fpoint, d, l)
    end
    function Track(ipoint::TPoint, d::Direction, b::Box)
        fpoint = intersect_box(ipoint, d, b)
        l = norm(fpoint-ipoint)
        new(ipoint, fpoint, d, l)
    end
end

function (t::Track)(λ)
    t.ipoint+λ*t.norm*t.direction
end

function intersect_box(p::TPoint, d::Direction, box::Box)
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

function Base.reverse(t::Track)
    Track(t.fpoint, t.ipoint)
end

function reduce_f(t::Track, f)
    g(λ) = f(t(λ).x, t(λ).y)
end

function inverse_column_depth(tr::Track, cd, valley)
    f(λ) = column_depth(tr, λ, valley) - cd
    λ_int = find_zero(f, (0,1))
    λ_int
end

function column_depth(t, Λ, valley; ρ_air=ρ_air0, ρ_rock=ρ_rock0)
#function column_depth(t::Track, Λ, valley; ρ_air=ρ_air0, ρ_rock=ρ_rock0)
    oned_valley = reduce_f(t, valley)
    root_func(λ) = oned_valley(λ)-t(λ).z
    zeros = find_zeros(root_func, 0, Λ)
    rgen = vcat([0], zeros, [Λ])
    ranges = [
        (x[2]-x[1], is_inside(t((x[1]+x[2])/2), valley))
        for x in zip(rgen[1:end-1], rgen[2:end])
    ]
    cd = 0mwe
    for x in ranges
        width, in_mountain = x
        ρ = in_mountain ? ρ_rock : ρ_air
        cd += width*t.norm*ρ
    end
    cd
end

#function total_column_depth(t::Track, valley; ρ_air=1.225e-3, ρ_rock=2.6)
function total_column_depth(t, valley; ρ_air=1.225e-3, ρ_rock=2.6)
    column_depth(t, 1, valley)
end

end # module
