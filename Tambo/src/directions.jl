struct Direction
    θ::Float64
    ϕ::Float64
    proj::SVector{3,Float64}
end

function Direction(θ::T, ϕ::U) where {T,U<:Number}
    θ, ϕ = Float64(θ), Float64(ϕ)
    proj = SVector{3}([cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)])
    return Direction(θ, ϕ, proj)
end

function Direction(x::T, y::U, z::V) where {T,U,V<:Number}
    x, y, z = promote(x, y, z)
    proj = SVector{3}([x, y, z]) ./ norm((x, y, z))
    θ = acos(proj.z)
    ϕ = atan(proj.y, proj.x)
    return Direction(θ, ϕ, proj)
end

function Direction(sv::SVector{3})
    return Direction(sv.x, sv.y, sv.z)
end

function Direction(ip::SVector{3}, fp::SVector{3})
    return Direction(ip .- fp)
end

function Direction(pp_vector::PyObject)
    return Direction(pp_vector.x, pp_vector.y, pp_vector.z)
end

function Base.show(io::IO, d::Direction)
    print(
        io,
        """
        θ (degrees): $(d.θ * 180 / π)°
        ϕ (degrees): $(d.ϕ * 180 / π)°
        proj: [$(d.proj.x), $(d.proj.y), $(d.proj.z)]
        """
    )
end

Base.reverse(d::Direction) = Direction(-d.proj...)

Base.:/(sv::SVector{3}, d::Direction) = sv ./ d.proj
Base.:*(m, d::Direction) = m .* d.proj

struct Box{T<:Number}
    c1::SVector{3,T}
    c2::SVector{3,T}
end