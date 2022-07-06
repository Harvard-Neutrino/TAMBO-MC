module ProposalInterface

using PyCall
using ProgressBars
using DelimitedFiles

export Particle, make_sector, define_particle, make_medium, make_propagator, propagate, vector3D

pp = pyimport("proposal")


mutable struct Particle
    pdg_mc::Int64
    Eᵢ::Float64
    position
    direction
    E_final
    range::Float64
    parent
    children::Array{Particle}
end

const propagators = Dict()

function make_sector(medium, start, stop)

    # Defining a Sector
    sec_def = pp.SectorDefinition()

    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), stop, start)

    if medium[1] == 1
        sec_def.medium = pp.medium.StandardRock(medium[2])
    else 
        sec_def.medium = pp.medium.Air([medium[2]])
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

function vector3D(vector)
    return pp.Vector3D(vector[1], vector[2], vector[3])
end

function make_medium(medium)

    # Definig sectors from medium (Array) and calculating detector length
    sectors = Vector()

    detector_length = 0.0
    sector_count = 0
    for (composition,l) in medium

        sector_count += 1

        # Calculate start and stop points for each sector
        start = detector_length
        # update detector length
        detector_length += l
        stop = detector_length

        push!(sectors,make_sector(composition,start,stop))
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
                        detector=pp.geometry.Sphere(pp.Vector3D(),detector_length, 0),
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
        elseif sec.type < 1000000001 
            # Should I worry about the angle? is there any way to get it?
            if sec.type in [13, 15, 11] 
                child = Particle(sec.type + 1, sec.energy, sec.position, sec.direction, 0.0, 0.0, parent_particle, [])
                push!(parent_particle.children, child)
            elseif sec.type in [-13, -15, -11] 
                child = Particle(sec.type - 1, sec.energy, sec.position, sec.direction, 0.0, 0.0, parent_particle, [])
                push!(parent_particle.children, child)
            end
        end

    end

    parent_particle.range = secondaries.particles[end].position.magnitude()

    E_final = Dict("continuous" => continuous, "epair" => epair, "brems" => brems, "ioniz" => ioniz, "photo" => photo)

    parent_particle.E_final = E_final

end


function propagate(particle::Particle, medium = nothing, propagator = nothing)

    particle_def = define_particle(particle)

    if propagator != nothing
        prop = propagator

    # Stores previous propagators to make code faster
    elseif particle.pdg_mc in keys(propagators)
        prop = propagators[particle.pdg_mc]
    else
        prop = make_propagator(particle,medium)
        propagators[particle.pdg_mc] = prop
    end

    lepton = pp.particle.DynamicData(particle_def.particle_type)
    # lepton.position = pp.Vector3D(0, 0, 0)
    # lepton.direction = pp.Vector3D(0, 0, -1)
    lepton.position = particle.position
    lepton.direction = particle.direction
    lepton.energy = particle.Eᵢ
    lepton.propagated_distance = 0.0
    lepton.time = 0.0

    secondaries = prop.propagate(lepton)

    analyze_secondaries!(secondaries, particle)

    if length(particle.children) > 0
        for child in particle.children
            propagate(child, medium)
        end
    end
end


function Base.show(io::IO, particle::Particle)
    
    print(io, 
    """{
        "Particle Type" : $(particle.pdg_mc),
        "Initial Energy" : $(particle.Eᵢ),
        "Range" : $(particle.range),
        "Energy breakdown" : $(particle.E_final),
        "Decay Products" : $(particle.children))
        }""")
end

end # module

##############################################################################################################
# EXAMPLE CODE

# using .ProposalInterface
# using Statistics
# using ProgressBars

# medium = [(0,1e7),(1,1e7)]
# energies = collect(6:0.5:11)

    
# for energy in energies
#     println("@ $energy MeV:")

#     for i in tqdm(1:1:200)

#         particle = Particle(14, 10.0^energy, vector3D([0 0 0]), vector3D([0 0 -1])  ,nothing, 0.0, nothing, [])
#         propagate(particle, medium)
#         # println(particle)
#     end
# end

##############################################################################################################
# MUON RANGES CODE

using .ProposalInterface
using Statistics
using ProgressBars
using DelimitedFiles

# medium = [((0,1.0),1e7),((1,1.0),1e7)]
medium = [((1,1.0),1e10)]
energies = collect(6:0.25:11)

results_dir = "/n/home08/jgarciaponce/Results/proposalinterface"
# results_file = "muon_range_stdrock.csv"



Results = Dict()
for energy in energies
    println("@ $energy MeV:")

    ranges = []

    for i in tqdm(1:1:1000)

        particle = Particle(14, 10.0^energy, vector3D([0 0 0]), vector3D([0 0 -1])  ,nothing, 0.0, nothing, [])
        propagate(particle, medium)
        push!(ranges, particle.range)

        # println(particle)
    end

    Results[energy] = ranges

    results_file = "muon_range_stdrock_@$(energy).csv"

    results_path = results_dir * "/" * results_file

    mean_range = mean(ranges)
    standard_deviation = std(ranges)

    writedlm(results_path,  ranges, ',')

    println("""

    Mean Range: $(mean_range)
    Sigma : $(standard_deviation)

    Saved to file: $(results_path)

    """)

end

