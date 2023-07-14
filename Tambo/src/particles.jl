struct Particle
    pdg_mc::Int64
    energy::Float64
    position::SVector{3, Float64}
    direction::Direction
    parent::Union{Particle,Nothing}
end

function Base.show(io::IO, particle::Particle)
    s = "pdg_mc=$(particle.pdg_mc), "
    s *= "energy=$(particle.energy / units.GeV) GeV, "
    s *= "position=$(particle.position / units.m)m, "
    s *= "direction=$(particle.direction)"
    print(io, s)
end

const range_parameters = Dict(
    14 => (1.76666667e-1 * units.GeV / units.mwe, 2.0916666667e-4 / units.mwe),
    16 => (1.473684210526e3 * units.GeV / units.mwe, 2.63e-5 / units.mwe),
    -14 => (1.76666667e-1 * units.GeV / units.mwe, 2.0916666667e-4 / units.mwe),
    -16 => (1.473684210526e3 * units.GeV / units.mwe, 2.63e-5 / units.mwe)
)

"""
    lepton_range(e, pdg_code)

Compute the 99.9% column depth for a lepton with energy `e` using parametrization
from https://doi.org/10.1016/j.cpc.2021.108018.

# Example
```julia-repl
julia> lepton_range(1units.PeV, false) / units.mwe

julia> lepton_range(1units.PeV, true) / units.mwe

```
"""
function lepton_range(e::Float64, pdg_code::Int)
    α, β = range_parameters[pdg_code]
    range = log(1 + e * β / α) / β
    return range
end

function Base.getindex(v::Vector{Particle}, s::String)
    return getfield.(v, Symbol(s))
end
