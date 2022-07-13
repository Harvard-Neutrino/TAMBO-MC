mutable struct Particle
    pdg_mc::Int64
    energy::Float64
    init_position::SVector{3}
    final_position::SVector{3}
    direction::SVector{2}
    losses::Vector{Vector}
    parent::Union{Particle, Nothing}
    children::Array{Particle}
end

"""
    lepton_range(e, is_tau)

Compute the 99.9% column depth for a lepton with energy `e` using parametrization
from https://doi.org/10.1016/j.cpc.2021.108018. If `is_tau` is true false this
is the range for a μ. If `is_tau` is true, the range for a τ is added to that of
a μ

# Example
```julia-repl
julia> lepton_range(1units.PeV, false) / units.mwe

julia> lepton_range(1units.PeV, true) / units.mwe

```
"""
function lepton_range(e::Float64, is_tau::Bool)   
    αμ = 1.76666667e-1 * units[:GeV] / units[:mwe]
    βμ = 2.0916666667e-4 / units[:mwe]
    range = log(1 + e * βμ / αμ) / βμ
    if is_tau
        ατ = 1.473684210526e3 * units[:GeV] / units[:mwe]
        βτ = 2.63e-5 / units[:mwe]
        range += log(1 + e * βτ / ατ) / βτ
    end
    range
end