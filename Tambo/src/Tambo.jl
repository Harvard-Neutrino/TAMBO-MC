__precompile__()
module Tambo

export Injector

using StaticArrays
using Rotations
using Random
using PyCall
using LinearAlgebra: norm
using Roots: find_zeros, find_zero
using Dierckx
using JLD2: jldopen
using ProgressBars
using Distributions: Uniform

const pp = PyNULL()
const TauMinusDef = PyNULL()
const TauPlusDef = PyNULL()
const MuMinusDef = PyNULL()
const MuPlusDef = PyNULL()
const EMinusDef = PyNULL()
const EPlusDef = PyNULL()
const EPlusAirCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const EMinusAirCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const MuPlusAirCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const MuMinusAirCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const TauPlusAirCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const TauMinusAirCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const EPlusRockCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const EMinusRockCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const MuPlusRockCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const MuMinusRockCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const TauPlusRockCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]
const TauMinusRockCross = [PyNULL(), PyNULL(), PyNULL(), PyNULL()]

function __init__()
    copy!(pp, pyimport("proposal"))
    pp.InterpolationSettings.tables_path = realpath(
        "$(@__DIR__)/../..//resources/proposal_tables/"
    )
    copy!(TauMinusDef, pp.particle.TauMinusDef())
    copy!(TauPlusDef, pp.particle.TauPlusDef())
    copy!(MuMinusDef, pp.particle.MuMinusDef())
    copy!(MuPlusDef, pp.particle.MuPlusDef())
    copy!(EMinusDef, pp.particle.EMinusDef())
    copy!(EPlusDef, pp.particle.EPlusDef())
    function make_pp_crosssection(particle_def, name)
        ecut = 50
        vcut = 1e-2
        cuts = pp.EnergyCutSettings(ecut, vcut, false)
        interpolate = true
        target = getproperty(pp.medium, name)()
        cross = pp.crosssection.make_std_crosssection(;
            particle_def=particle_def, target=target, interpolate=interpolate, cuts=cuts
        )

        return cross
    end
    xs = make_pp_crosssection(EMinusDef, "Air")
    ys = EMinusAirCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(EPlusDef, "Air")
    ys = EPlusAirCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(MuMinusDef, "Air")
    ys = MuMinusAirCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(MuPlusDef, "Air")
    ys = MuPlusAirCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(TauMinusDef, "Air")
    ys = TauMinusAirCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(TauPlusDef, "Air")
    ys = TauPlusAirCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(EMinusDef, "StandardRock")
    ys = EMinusRockCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(EPlusDef, "StandardRock")
    ys = EPlusRockCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(MuMinusDef, "StandardRock")
    ys = MuMinusRockCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(MuPlusDef, "StandardRock")
    ys = MuPlusRockCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(TauMinusDef, "StandardRock")
    ys = TauMinusRockCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
    xs = make_pp_crosssection(TauPlusDef, "StandardRock")
    ys = TauPlusRockCross
    for (x, y) in zip(xs, ys)
        copy!(y, x)
    end
end

include("units.jl")
include("locations.jl")
include("directions.jl")
include("geometries.jl")
include("tracks.jl")
#include("segments.jl")
include("particles.jl")
include("inject.jl")
include("proposal.jl")

@Base.kwdef mutable struct Simulator
    n::Int = 10
    geo_spline_path::String = realpath("$(@__DIR__)/../../resources/tambo_spline.jld2")
    diff_xs_path::String = realpath(
        "$(@__DIR__)/../../resources/cross_sections/tables/csms_differential_cdfs.h5"
    )
    tambo_coordinates::Coord = minesite_coord
    ν_pdg::Int = 16
    γ::Float64 = 1
    emin::Float64 = 1e6units.GeV
    emax::Float64 = 1e9units.GeV
    θmin::Float64 = 0.0
    θmax::Float64 = π
    ϕmin::Float64 = 0.0
    ϕmax::Float64 = 2π
    r_injection::Float64 = 900units.m
    l_endcap::Float64 = 1units.km
    seed::Int64 = 0
end

end # module
