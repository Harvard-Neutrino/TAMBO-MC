struct UniformAngularSampler
    θmin::Float64
    θmax::Float64
    ϕmin::Float64
    ϕmax::Float64
    function UniformAngularSampler(θmin, θmax, ϕmin, ϕmax)
        ϕmin = mod(ϕmin, 2π)
        ϕmax = mod(ϕmax, 2π)
        @assert ϕmin <= ϕmax "ϕmin greater than ϕmax"
        @assert θmin <= θmax "θmin greater than θmax"
        return new(θmin, θmax, ϕmin, ϕmax)
    end
end

function Base.show(io::IO, uniformsampler::UniformAngularSampler)
    s = "UniformAngularSampler("
    s *= "θmin=$(uniformsampler.θmin / π * 180)∘, "
    s *= "θmax=$(uniformsampler.θmax / π * 180)∘, "
    s *= "ϕmin=$(uniformsampler.ϕmin / π * 180)∘, "
    s *= "ϕmax=$(uniformsampler.ϕmax / π * 180)∘"
    s *= ")"
    print(io, s)
end

"""
    Base.rand(uniformsampler::UniformAngularSampler)

Return zenith (θ) and azimuth (ϕ) angles sampled uniformly on a sphere

TBW
"""
function Base.rand(uniformsampler::UniformAngularSampler)
    # Randomly sample zenith uniform in phase space
    θ = acos(rand(Uniform(cos(uniformsampler.θmax), cos(uniformsampler.θmin))))
    # Randomly sample azimuth
    ϕ = rand(Uniform(uniformsampler.ϕmin, uniformsampler.ϕmax))
    return θ, ϕ
end

function probability(sampler::UniformAngularSampler)
    Ω = (sampler.ϕmax - sampler.ϕmin) * (cos(sampler.θmin) - cos(sampler.θmax))
    return 1 / Ω
end