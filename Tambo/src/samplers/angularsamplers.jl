struct UniformAngularSampler
    θmin::Float64
    θmax::Float64
    ϕmin::Float64
    ϕmax::Float64
end

function Base.show(io::IO, uniformsampler::UniformAngularSampler)
    print(
        io,
        """
        θmin (degrees): $(uniformsampler.θmin / π * 180)
        θmax (degrees): $(uniformsampler.θmax / π * 180)
        ϕmin (degrees): $(uniformsampler.ϕmin / π * 180)
        ϕmax (degrees): $(uniformsampler.ϕmax / π * 180)
        """,
    )
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
