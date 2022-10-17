using PyCall
using StaticArrays

"""
This function wrapping is necessary to prevent PROPOSAL being
a PyNull object. I don't know what that means or why this stops it
https://discourse.jullialang.org/t/can-we-import-dill-in-pycall-rather-than-pickle
As per this discussion, pyimports are cached so this shouldn't cause
a performance hit
"""

const particle_def_dict = Dict(
    11 => EMinusDef,
    -11 => EPlusDef,
    13 => MuMinusDef,
    -13 => MuPlusDef,
    15 => TauMinusDef,
    -15 => TauPlusDef,
)

const pp_crosssections_dict = Dict(
    (11, "Air") => EMinusAirCross,
    (-11, "Air") => EPlusAirCross,
    (13, "Air") => MuMinusAirCross,
    (-13, "Air") => MuPlusAirCross,
    (15, "Air") => TauMinusAirCross,
    (-15, "Air") => TauPlusAirCross,
    (11, "StandardRock") => EMinusRockCross,
    (-11, "StandardRock") => EPlusRockCross,
    (13, "StandardRock") => MuMinusRockCross,
    (-13, "StandardRock") => MuPlusRockCross,
    (15, "StandardRock") => TauMinusRockCross,
    (-15, "StandardRock") => TauPlusRockCross,
)

struct Loss
    int_type::Int64
    energy::Float64
    position::SVector{3,Float64}
    function Loss(i, e, p)
        return new(i, e, p)
    end
end

function Loss(int_type::Int, e::Float64, position::PyObject)
    position = SVector{3}([position.x, position.y, position.z])
    return Loss(int_type, e, position)
end

function Base.show(io::IO, l::Loss)
    return print(
        io,
        """
        "Interaction Type" : $(l.int_type),
        "Energy (GeV)" : $(l.energy / units.GeV)
        "Position (m)" : $(l.position ./ units.m)
        """,
    )
end

struct ProposalResult
    losses::Vector{Loss}
    did_decay::Bool
    decay_products::Vector{Particle}
    final_pos::SVector{3}
    final_lepton_energy::Float64
end

function ProposalResult(pos)
    losses = Loss[]
    did_decay = false
    decay_products = Particle[]
    return ProposalResult(losses, did_decay, decay_products, pos)
end

function ProposalResult(secondaries, parent_particle)
    losses = Loss[]
    for sec in secondaries.stochastic_losses()
        int_type = sec.type
        sec_e = sec.energy
        sec_e = sec.energy * units.MeV
        pos = position_from_pp_vector(sec.position)
        l = Loss(int_type, sec_e, pos)
        push!(losses, l)
    end

    children = Particle[]
    for product in secondaries.decay_products()
        pdg_code = product.type
        energy = product.energy * units.MeV
        position = position_from_pp_vector(product.position)
        direction = Direction(product.direction)
        parent = parent_particle
        child = Particle(pdg_code, energy, position, direction, parent)
        push!(children, child)
    end
    did_decay = length(children) > 0
    final_pos = position_from_pp_vector(secondaries.final_state().position)
    final_e = secondaries.final_state().energy * units.MeV
    return ProposalResult(losses, did_decay, children, final_pos, final_e)
end

function show(io::IO, result::ProposalResult)
    return print(
        io,
        """
        losses:
        $(result.losses)

        decay_products: 
        $(result.decay_products)

        did_decay: $(result.did_decay)
        final_pos (m): $(result.final_pos / units.m)
        final_lepton_energy (GeV): $(result.final_lepton_energy / units.GeV)
        """,
    )
end

function pp_particle_def(particle::Particle)
    return particle_def_dict[particle.pdg_mc]
end

function position_from_pp_vector(pp_vector)
    return SVector{3}([pp_vector.x, pp_vector.y, pp_vector.z]) .* units.cm
end

"""
    make_pp_vector(v::SVector{3})

TBW
"""
function make_pp_vector(v::SVector{3})
    v = v ./ units.cm
    return pp.Cartesian3D(v)
end

function make_pp_direction(d::Direction)
    return pp.Cartesian3D(d.proj)
end

function make_pp_crosssection(particle::Particle, medium_name)
    return pp_crosssections_dict[(particle.pdg_mc, medium_name)]
end

function make_pp_utility(particle_def, cross)
    collection = pp.PropagationUtilityCollection()

    collection.displacement = pp.make_displacement(cross, true)
    collection.interaction = pp.make_interaction(cross, true)
    collection.time = pp.make_time(cross, particle_def, true)
    collection.decay = pp.make_decay(cross, particle_def, true)

    utility = pp.PropagationUtility(; collection=collection)
    return utility
end

function make_pp_geometry(start, end_)
    return pp.geometry.Sphere(pp.Cartesian3D(), end_ / units.cm, start / units.cm)
end

function make_pp_density_distribution(density)
    return pp.density_distribution.density_homogeneous(density / (units.gr / units.cm^3))
end

function make_propagator(
    particle::Particle,
    particle_def,
    media::Vector{String},
    densities::Vector{Float64},
    lengths::Vector{Float64},
)
    geometries = []
    start = 0
    push!(lengths, 1e10 * units.km)
    for l in lengths
        end_ = start + l
        geo = make_pp_geometry(start, end_)
        push!(geometries, geo)
        start = end_
    end
    density_distributions = [make_pp_density_distribution(d) for d in densities]
    push!(
        density_distributions,
        make_pp_density_distribution(units.ρair0 / (units.gr / units.cm^3)),
    )
    crosses = [make_pp_crosssection(particle, m) for m in media]
    push!(crosses, make_pp_crosssection(particle, "Air"))
    utilities = [make_pp_utility(particle_def, cross) for cross in crosses]
    dumb_list = [x for x in zip(geometries, utilities, density_distributions)]
    prop = pp.Propagator(particle_def, dumb_list)
    return prop
end

function propagate_charged_lepton(
    clepton::Particle,
    media::Vector{String},
    densities::Vector{Float64},
    lengths::Vector{Float64},
)
    particle_def = pp_particle_def(clepton)

    prop = make_propagator(clepton, particle_def, media, densities, lengths)

    lepton = pp.particle.ParticleState()
    lepton.position = make_pp_vector(clepton.position)
    lepton.direction = make_pp_direction(clepton.direction)
    lepton.energy = clepton.energy / units.MeV
    lepton.propagated_distance = 0.0
    lepton.time = 0.0

    secondaries = prop.propagate(lepton)
    return secondaries
end

function propagate(v::Vector{Particle}, geo::Geometry; track_progress=true)
    if track_progress
        iter = ProgressBar(v)
    else
        iter = v
    end
    return [propagate(p, geo) for p in iter]
end

function propagate(final_state::Particle, geo::Geometry)
    # Reverse Direction since Track tells us where we're going
    # But Particle.direction tells us where it is from
    t = Track(final_state.position, reverse(final_state.direction), geo.box)
    ranges = computeranges(t, geo)
    # I feel like ranges is an artificial object... not gonna fix now though
    # TODO refactor ranges
    result = propagate(
        final_state,
        getfield.(ranges, :medium_name),
        getfield.(ranges, :density),
        getfield.(ranges, :length),
    )
    return result
end

function propagate(
    particle::Particle,
    media::Vector{String},
    densities::Vector{Float64},
    lengths::Vector{Float64},
)
    if abs(particle.pdg_mc) in [11, 13, 15]
        secondaries = propagate_charged_lepton(particle, media, densities, lengths)
        result = ProposalResult(secondaries, particle)
        # We have decided to let CORSIKA do the handling of decay products
        # This will need a little rethinking if we change on that
        #for child in result.decay_products
        #    propagate(child, medium)
        #end
    else
        result = ProposalResult(particle.position)
    end
    return result
end

function propagate(
    vp::Vector{Particle},
    media::Vector{String},
    densities::Vector{Float64},
    lengths::Vector{Float64};
    track_progress=true,
)
    if track_progress
        iter = ProgressBar(vp)
    else
        iter = vp
    end
    return [propagate(p, media, densities, lengths) for p in iter]
end
