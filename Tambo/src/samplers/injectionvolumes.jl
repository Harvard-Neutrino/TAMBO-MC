using StaticArrays: SVector

struct SymmetricInjectionCylinder
    r_injection::Float64
    l_endcap::Float64
end

function Base.show(io::IO, cylinder::SymmetricInjectionCylinder)
    print(
        io,
        """
        r_injection (m): $(cylinder.r_injection / units.m)
        l_endcap (m): $(cylinder.l_endcap / units.m)
        """
    )
end

struct AsymmetricInjectionCylinder
    r_injection::Float64
    l_incoming_endcap::Float64
    l_outgoing_endcap::Float64
end

function Base.show(io::IO, cylinder::AsymmetricInjectionCylinder)
    print(
        io,
        """
        r_injection (m): $(cylinder.r_injection / units.m)
        l_incoming_endcap (m): $(cylinder.l_incoming_endcap / units.m)
        l_outgoing_endcap (m): $(cylinder.l_outgoing_endcap / units.m)
        """
    )
end

"""
    Base.rand(cylinder::SymmetricInjectionCylinder)

Sample an point of closest approach uniformly on a circle lying
in the xy-plane

# Example
```julia-repl
julia> injector = SymmetricInjectionCylinder(900units.m, 1units.km);
julia> rand.seed!(0);
julia> rand(injector)
TBW
```
"""
function Base.rand(cylinder::SymmetricInjectionCylinder)
    b = cylinder.r_injection .* sqrt(rand())
    ψ = rand(Uniform(0, 2π))
    p = SVector{3}([b * cos(ψ), b * sin(ψ), 0])
    return p
end

"""
    Base.rand(cylinder::AsymmetricInjectionCylinder)

Sample an point of closest approach uniformly on a circle of
radius r = cylinder.r_injection lying in the xy-plane. This can then
be rotated to a plane perpendicular to a particle's direction to
give the point of closest approach

# Example
```julia-repl
julia> injector = AsymmetricInjectionCylinder(900units.m, 1units.km);
julia> rand.seed!(0);
julia> rand(injector)
TBW
```
"""
function Base.rand(cylinder::AsymmetricInjectionCylinder)
    b = cylinder.r_injection .* sqrt(rand())
    ψ = rand(Uniform(0, 2π))
    p = SVector{3}([b * cos(ψ), b * sin(ψ), 0])
    return p
end

"""
    (cylinder::Injectioncylinder)(p::SVector{3})

TBW
"""
function probability(cylinder::SymmetricInjectionCylinder)
    return 1 / (π * cylinder.r_injection^2)
end