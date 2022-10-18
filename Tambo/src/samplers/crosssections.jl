using HDF5: h5open
using Dierckx: Spline2D

include("../units.jl")

const pdg_to_str = Dict{Int,String}(
    12 => "nue",
    14 => "numu",
    16 => "nutau",
    -12 => "nuebar",
    -14 => "numubar",
    -16 => "nutaubar",
)

struct OutgoingCCEnergy
    fname::String
    spline::Spline2D
    _emin_z::Float64
    _es::Vector{Float64}
end

function OutgoingCCEnergy(fname::String, ν_pdg::Int)
    nutype = ν_pdg > 0 ? "nu" : "nubar"
    h5f = h5open(fname)
    log10_es = h5f["log10_es"][:]
    # This number is hardcoded in how we generated the cross cross_sections
    # If we start generating XS with anything besides
    # z = (eout - 10 GeV) / (ein - 10 GeV) we need to change this
    _emin_z = 10.0 * units.GeV
    spl = Spline2D(log10_es, h5f["cdf_vals"][:], h5f["$(nutype)_cc_cdf"][:, :]; kx=1, ky=1)
    close(h5f)
    return OutgoingCCEnergy(realpath(fname), spl, _emin_z, 10.0 .^ (log10_es))
end

function Base.show(io::IO, sampler::OutgoingCCEnergy)
    print(
        io,
        """
        fname: $(sampler.fname)
        emin (GeV): $(sampler._es[1] / units.GeV)
        emax (GeV): $(sampler._es[end] / units.GeV)
        """,
    )
end

function Base.rand(sampler::OutgoingCCEnergy, ein::Float64)
    z = 1 - sampler.spline(log(10.0, ein), rand())
    eout = z * (ein - sampler._emin_z) + sampler._emin_z
    return eout
end

function Base.rand(n::Int, sampler::OutgoingCCEnergy, ein::Float64)
    return [rand(sampler, ein) for _ in 1:n]
end
