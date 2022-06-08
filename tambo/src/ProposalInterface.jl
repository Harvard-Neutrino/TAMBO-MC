module ProposalInterface

using PyCall
using ProgressBars
using DelimitedFiles

export Particle, make_sector, define_particle, make_medium, make_propagator, propagate

pp = pyimport("proposal")


# Make this useful
mutable struct Particle
    pdg_mc::Int64
    energy::Float64
    E
    range::Float64
    parent
    children::Array{Particle}
end

function make_sector(medium, start, stop)

    # Defining a Sector
    sec_def = pp.SectorDefinition()

    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), stop, start)

    if medium == 1
        sec_def.medium = pp.medium.StandardRock()

    else 
        sec_def.medium = pp.medium.Air()
    end

    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere

    sec_def.crosssection_defs.brems_def.lpm_effect = true
    sec_def.crosssection_defs.epair_def.lpm_effect = true
    sec_def.cut_settings.ecut = -1.0
    sec_def.cut_settings.vcut = 1e-3

    sec_def.do_continuous_randomization = true
    sec_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.AbramowiczLevinLevyMaor97

    return sec_def

end


function define_particle(particle::Particle)
    particles = Dict(
        16 => :(pp.particle.TauMinusDef()),
        -16 => :(pp.particle.TauPlusDef()),
        14 => :(pp.particle.MuMinusDef()),
        -14 => :(pp.particle.MuPlusDef()),
        12 => :(pp.particle.EMinusDef()),
        -12 => :(pp.particle.EPlusDef())
    )
    # {16: 'TauMinusDef', -16: 'TauPlusDef',
    #          14: 'MuMinusDef',  -14: 'MuPlusDef',
    #          12: 'EMinusDef',   -12: 'EPlusDef',}

    return eval(particles[particle.pdg_mc])
end


function make_medium(medium)

    # Definig sectors from medium (Array) and calculating detector length
    sectors = Vector()

    # println("Defining Sectors:")

    detector_length = 0.0
    sector_count = 0
    for (rock,l) in medium

        sector_count += 1

        # Calculating start and stop points for each sector
        start = detector_length
        # updates detector length
        detector_length += l
        stop = detector_length

        push!(sectors,make_sector(rock,start,stop))
    end

    return sectors, detector_length
end

function make_propagator(particle, medium)

    sectors, detector_length = make_medium(medium)

    particle_def = define_particle(particle)


    # Interpolation tables
    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "/tmp/tables"
    interpolation_def.path_to_tables_readonly = "/tmp/tables"


    propagator = pp.Propagator(particle_def=particle_def, 
                        sector_defs=sectors,
                        detector=pp.geometry.Sphere(pp.Vector3D(),detector_length, 0),# How do I choose the radius? Sum of Lengths?
                        interpolation_def=interpolation_def)
    return propagator
    
end

function analyze_secondaries!(secondaries, parent_particle)

    continuous = Vector{Vector}()
    epair = Vector{Vector}()
    brems = Vector{Vector}()
    ioniz = Vector{Vector}()
    photo = Vector{Vector}()

    for sec in secondaries.particles

        log_sec_energy = log.(10,sec.parent_particle_energy .- sec.energy)
        
        if sec.type == 1000000001
            push!(continuous, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
        elseif sec.type == 1000000004
            push!(epair, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
        elseif sec.type == 1000000002
            push!(brems, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
        elseif sec.type == 1000000003
            push!(ioniz, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
        elseif sec.type == 1000000005
            push!(photo,[log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
        elseif sec.type < 1000000001 && sec.type in [13, 15, -13, -15]
            # Should I worry about the angle? is there any way to get it?
            child = Particle(sec.type, sec.energy, 0.0, 0.0, parent_particle, [])
            push!(parent_particle.children, child)
        end

    end

    parent_particle.range = secondaries.particles[end].position.magnitude()

    E = Dict("continuous" => continuous, "epair" => epair, "brems" => brems, "ioniz" => ioniz, "photo" => photo)

    parent_particle.E = E
    # println(parent_particle.E)
    # readline()

end


function propagate(particle::Particle, medium = nothing, propagator = nothing)

    particle_def = define_particle(particle)

    if propagator == nothing
        prop = make_propagator(particle,medium)
    else
        prop = propagator
    end

    # Probably gonna have to fix this.
    lepton = pp.particle.DynamicData(particle_def.particle_type)
    lepton.position = pp.Vector3D(0, 0, 0)
    lepton.direction = pp.Vector3D(0, 0, -1)
    lepton.energy = particle.energy
    # lepton.energy = 1e7
    lepton.propagated_distance = 0.0
    lepton.time = 0.0

    println("\n\nSTARTING PROPAGATION\n\n")

    secondaries = prop.propagate(lepton)

    analyze_secondaries!(secondaries, particle)

    if length(particle.children) > 0
        for child in particle.children
            propagate(child, medium)

        end
    end

end


function Base.show(io::IO, particle::Particle)
    # base case
    
    print(io, 
    """{
        "Particle Type" : $(particle.pdg_mc),
        "Initial Energy" : $(particle.energy),
        "Range" : $(particle.range),
        "Energy breakdown" : $(particle.E),
        "Decay Products" : $(particle.children))
        }""")
end

end # module


























##############################################################################################################
# Muon Range Code

using .ProposalInterface
using Statistics
using DelimitedFiles
using JSON

particle = Particle(12, 10.0^8, nothing, 0.0, nothing, [Particle(14, 10.0^8, nothing, 0.0, nothing, [Particle(14, 10.0^8, nothing, 0.0, nothing, [])])])

particle.children[1].parent = particle
particle.children[1].children[1].parent = particle.children[1]


# test
medium = [(0,1e7),(1,1e7)]
# energy = collect(6:0.5:11)

propagate(particle, medium)

display(particle)




# function recursive_print(particle::Particle)

#     if particle.children == []
#         event = Dict("type" => particle.pdg_mc, "range" => particle.range, "E" => particle.E)
#         println(event)

#         if particle.parent != nothing
#             particle.parent.children = setdiff(particle.parent.children, [particle])
#             recursive_print(particle.parent)
#         end

#     else
#         # broadcast(recursive_print, particle.children)
#         # recursive_print.(particle.children)
#         for child in particle.children
#             recursive_print(child)
#         end
#     end
# end

# recursive_print(particle)

# energies = Vector{Float64}()
# range = Vector{Float64}()
# error = Vector{Float64}()

# for e in energy
#     println("\n\n\nLogEnergy: $e MeV")
#     lepton_length, n_secondaries, E = ProposalInterface.propagate_mcp(particle, 10^e, 1000, medium, prop)
#     d = mean(lepton_length)
#     # d_logs = log.(10,lepton_length./1000)
#     # d = mean(d_logs)
#     err = std(lepton_length)
#     # err = std(d_logs)
#     println("Mean Range: $d")
#     println("Standard Deviation: $err")
#     push!(range, d)
#     push!(energies, e)
#     push!(error, err)
# end

# println(error)
# println(range)
# println(energies)

# writedlm( "/n/home08/jgarciaponce/Results/Range.csv",  range, ',')
# writedlm( "/n/home08/jgarciaponce/Results/Error.csv",  error, ',')
# writedlm( "/n/home08/jgarciaponce/Results/Energies.csv",  energies, ',')
# writedlm( "Results/leptonMinus_Secondaries.csv",  n_secondaries, ',')

# println("Files saved in Results Folder")

# lepton_length, n_secondaries, E = propagate_mcp(particle, medium, mcp_energies)

# print(lepton_length)
# print(n_secondaries)
# print(E)










# for x in 1:dims
#     particulas[x] = Particle(14,1.0, nothing, nothing, [])
# end
