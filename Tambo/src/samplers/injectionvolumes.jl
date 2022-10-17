using StaticArrays: SVector

struct InjectionVolume
    r_injection::Float64
    l_endcap::Float64
end

function Base.show(io::IO, volume::InjectionVolume)
    return print(
        io,
        """
        r_injection (m): $(volume.r_injection / units.m)
        l_endcap (m): $(volume.l_endcap / units.m)
        """,
    )
end

"""
    Base.rand(volume::InjectionVolume)

Sample an point of closest approach uniformly on a circle lying
in the xy-plane

# Example
```julia-repl
julia> injector = InjectionVolume(900units.m, 1units.km);
julia> rand.seed!(0);
julia> rand(injector)
TBW
```
"""
function Base.rand(volume::InjectionVolume)
    # Sample impact parameter uniformly on a disc
    b = volume.r_injection .* sqrt(rand())
    # Sample angle on disc 
    ψ = rand(Uniform(0, 2π))
    p = SVector{3}([b * cos(ψ), b * sin(ψ), 0])
    return p
end

"""
    (volume::InjectionVolume)(p::SVector{3})

TBW
"""
function (volume::InjectionVolume)(p::SVector{3})
    return 1 / (π * volume.r_injection^2)
end
