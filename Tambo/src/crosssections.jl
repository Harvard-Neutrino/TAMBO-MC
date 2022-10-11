using HDF5
using StatsBase
using GridInterpolations

const pdg_to_str = Dict{Int,String}(12 => "nue", 14 => "numu", 16 => "nutau")

struct DifferentialXS
    # Discuss if we want to allow for there not to be only one type of nucleon
    zs::LinRange
    a::Vector{Float64}
    es::Vector{Float64}
    emin::Float64
    xs::Matrix{Float64}
    function DifferentialXS(zs, a, es, xs)
        size(xs)[1] == size(es)[1] || error("We need a cross section for each energy")
        size(xs)[2] == size(zs)[1] || error("Number of zs wrong")
        size(xs)[2] == size(a)[1] || error("holder array is not the right size")
        zs[1] == 0 && zs[end] == 1 || error("zs need to go from 0 to 1")
        zs = LinRange(0, 1, length(a))
        emin = 1e10
        return new(zs, a, es, emin, xs)
    end
end

function DifferentialXS(fname::String, ν_pdg::Int, interaction::String)
    h5f = h5open(fname)
    es = 10 .^ h5f["energies"][:]
    zs = h5f["zs"][:]
    flavor = pdg_to_str[abs(ν_pdg)]
    nt = ν_pdg > 0 ? "nu" : "nubar"
    xs = h5f["$(flavor)_$(nt)_$(interaction)"][:, :]
    a = zeros(length(zs))
    close(h5f)
    return DifferentialXS(zs, a, es, xs)
end

function DifferentialXS()
    fname = "$(@__DIR__)/../../../resources/cross_sections/tables/csms_differential.h5"
    ν_pdg = 16
    interaction = "CC"
    return DifferentialXS(fname, ν_pdg, interaction)
end

function DifferentialXS(fname::String)
    ν_pdg = 16
    interaction = "CC"
    return DifferentialXS(fname, ν_pdg, interaction)
end

function DifferentialXS(fname::String, interaction::String)
    ν_pdg = 16
    return DifferentialXS(fname, ν_pdg, interaction)
end

function DifferentialXS(fname::String, ν_pdg::Int)
    interaction = "CC"
    return DifferentialXS(fname, ν_pdg, interaction)
end

struct TotalXS end

"""
    searchsortednearest(a, x)

Find the index element of sorted array `a` whose value is closest to `x`. If `x`
is equally close to two elements, the first index will be taken.

```julia-repl
julia> searchsortednearest([1, 3, 5, 10], 5)
3

julia> searchsortednearest([1, 3, 5, 10], 4)
2

```
"""
function searchsortednearest(a, x)
    idx = searchsortedfirst(a, x)
    if (idx == 1)
        return idx
    end
    if (idx > length(a))
        return length(a)
    end
    if (a[idx] == x)
        return idx
    end
    if (abs(a[idx] - x) < abs(a[idx - 1] - x))
        return idx
    else
        return idx - 1
    end
end

"""
    get_xs(a, e, es, xs)
Fills input array `a` with the the row of `xs` corresponding to where `e` is 
closest to `es`. 

This filling is actually faster than just grabbing the row directly.
"""
function get_xs!(a, e, es, xs)
    idx = searchsortednearest(es, e)
    for i in 1:length(a)
        a[i] = xs[idx, i]
    end
end

"""
    StatsBase.sample(n, xs, e)

Draw samples from 
"""
function StatsBase.sample(n::Int, xs::DifferentialXS, e::T) where {T<:Number}
    get_xs!(xs.a, e, xs.es, xs.xs)
    x = sample(xs.zs, Weights(xs.a), n)
    return x .* (e - xs.emin) .+ xs.emin
end

function StatsBase.sample(xs::DifferentialXS, e::T) where {T<:Number}
    return sample(1, xs, e)[1]
end
