using HDF5: h5open
using Interpolations: LinearInterpolation, Extrapolation

const pdg_to_str = Dict{Int,String}(
    12 => "nue_nu",
    14 => "numu_nu",
    16 => "nutau_nu",
    -12 => "nue_nubar",
    -14 => "numu_nubar",
    -16 => "nutau_nubar",
)

@enum Interaction begin
    CC = 1
    NC = 2
end

struct TotalXS
    fname::String
    interaction::Interaction
    itp::Extrapolation
end

struct DifferentialXS
    fname::String
    interaction::Interaction
    itp::Extrapolation
    log_emin_z::Float64
end

struct OutgoingEnergy
    fname::String
    interaction::Interaction
    itp::Extrapolation
    log_emin_z::Float64
end

struct CrossSection
    ν_pdg::Int
    interaction::Interaction
    energy_sampler::OutgoingEnergy
    total_xs::TotalXS
    differential_xs::DifferentialXS
end

function CrossSection(dir::String, model_name::String, ν_pdg::Int, interaction::Interaction)
    
    sampling_file = "$(dir)/$(model_name)_differential_cross_section_cdf.h5"
    total_file = "$(dir)/$(model_name)_total_cross_section.h5"
    differential_file = "$(dir)/$(model_name)_differential_cross_section.h5"
    tot = TotalXS(total_file, ν_pdg, interaction)
    diff = DifferentialXS(differential_file, ν_pdg, interaction)
    sampl = OutgoingEnergy(sampling_file, ν_pdg, interaction)
    return CrossSection(ν_pdg, interaction, sampl, tot, diff)
end

function TotalXS(fname::String, ν_pdg::Int, interaction::Interaction)
    nutype = pdg_to_str[ν_pdg]
    int_str = "CC"
    if interaction==Interaction(2)
        int_str = "NC"
    end
    h5f = h5open(fname)
    log_xs = h5f["log_$(nutype)_$(int_str)"][:]
    log_energys = h5f["log_energies"][:]
    close(h5f)
    itp = LinearInterpolation(log_energys, log_xs)
    return TotalXS(fname, interaction, itp)
end

function DifferentialXS(fname::String, ν_pdg::Int, interaction::Interaction)
    nutype = pdg_to_str[ν_pdg]
    int_str = "CC"
    if interaction==Interaction(2)
        int_str = "NC"
    end
    h5f = h5open(fname)
    log_xs = h5f["log_$(nutype)_$(int_str)"][:, :]
    log_energys = h5f["log_energies"][:]
    log_emin = read(h5f["log_energymin"])
    zs = h5f["zs"][:]
    close(h5f)
    itp = LinearInterpolation((log_energys, zs), log_xs)
    return DifferentialXS(fname, interaction, itp, log_emin)
end

function OutgoingEnergy(fname::String, ν_pdg::Int, interaction::Interaction)
    nutype = pdg_to_str[ν_pdg]
    int_str = "CC"
    if interaction==Interaction(2)
        int_str = "NC"
    end
    h5f = h5open(fname)
    log_energys = h5f["log_energies"][:]
    log_emin = read(h5f["log_energymin"])
    vs = h5f["$(nutype)_$(int_str)"][:, :]
    us = h5f["us"][:]
    close(h5f)
    itp = LinearInterpolation((log_energys, us), vs)
    return OutgoingEnergy(fname, interaction, itp, log_emin)
end

function (xs::TotalXS)(e)
    # TODO fix this hack
    if e <= 105 * units.GeV
        return 0
    end
    return 10 ^ xs.itp(log(10, e))
end

function (xs::DifferentialXS)(ein, eout)
    emin = 10^xs.log_emin_z
    if eout < emin
        return 0.0
    end
    z = (eout - emin) / (ein - emin)
    logein = log(10, ein)
    return 10 ^ xs.itp(logein, z)
end

function (xs::OutgoingEnergy)(ein, u)
    return xs.itp(log(10, ein), u)
end

function Base.rand(xs::OutgoingEnergy, ein::Float64)
    # Hack to catch love energy things that are beneath our tables
    if log(10, ein) < 11.000000000000433
        return 0.0
    end
    u = rand()
    eout = xs(ein, u) * (ein - 10^xs.log_emin_z) + 10^xs.log_emin_z
    return eout
end

function Base.rand(xs::CrossSection, ein::Float64)
    return rand(xs.energy_sampler, ein)
end

function Base.rand(n::Int, xs::OutgoingEnergy, ein::Float64)
    return [rand(xs, ein) for _ in 1:n]
end

function Base.show(io::IO, xs::TotalXS)
    s = "TotalXS(fname=$(xs.fname), interaction=$(xs.interaction))"
    print(io, s)
end

function Base.show(io::IO, xs::DifferentialXS)
    s = "DifferentialXS(fname=$(xs.fname), interaction=$(xs.interaction), emin=$(10^xs.log_emin_z / 1e9) GeV)"
    print(io, s)
end

function Base.show(io::IO, xs::OutgoingEnergy)
    s = "OutgoingEnergy(fname=$(xs.fname), interaction=$(xs.interaction), emin=$(10^xs.log_emin_z / 1e9) GeV)"
    print(io, s)
end

function probability(xs::CrossSection, ein::Float64, eout::Float64)
    @assert eout <= ein "Energy out of range"
    σ = xs.total_xs(ein)
    diff_σ = xs.differential_xs(ein, eout)
    return diff_σ / σ
end

function probability(xs::CrossSection, event::InjectionEvent)
    return probability(xs, event.entry_state.energy, event.final_state.energy)
end
