module ProposalInterface

using PyCall
using ProgressBars
using DelimitedFiles

export Particle, make_sector, define_particle, make_medium, make_propagator, propagate

pp = pyimport("proposal")


# Make this useful
mutable struct Particle
    pdg_mc::Int64
    energy
    E
    parent
    children
end

function make_sector(medium, start, stop)

    # Defining a Sector
    sec_def = pp.SectorDefinition()

    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), stop, start)

    if medium == 1
        sec_def.medium = pp.medium.StandardRock()
        # leptonon Test only
        # sec_def.medium = pp.medium.Ice()
        # println("Ice")
    else 
        sec_def.medium = pp.medium.Air()
    end

    # What does this mean? Is this right?
    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere

    # The way it was
    sec_def.crosssection_defs.brems_def.lpm_effect = true
    sec_def.crosssection_defs.epair_def.lpm_effect = true
    sec_def.cut_settings.ecut = -1.0
    sec_def.cut_settings.vcut = 1e-3
    # sec_def.crosssection_defs.brems_def.lpm_effect = false
    # sec_def.crosssection_defs.epair_def.lpm_effect = false
    # sec_def.cut_settings.ecut = 500
    # sec_def.cut_settings.vcut = 0.05

    # Again, not sure what this means
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

    println("Defining Sectors:")

    detector_length = 0.0
    sector_count = 0
    for (rock,l) in medium

        sector_count += 1
        # detector_length
        println("Sector $sector_count")
        println("Rock? $rock")
        println("Length: $l")

        # Calculating start and stop points for each sector
        start = detector_length
        # updates detector length
        detector_length += l
        stop = detector_length
        print(medium)
        push!(sectors,make_sector(rock,start,stop))
    end

    #println(sectors)
    println("Detector Length: $detector_length")

    return sectors, detector_length
end

function make_propagator(particle, medium)

    sectors, detector_length = make_medium(medium)

    # particle = Particle(1.0,1.0)

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

function analyze_secondaries(secondaries, parent)
    # decay_products = [p for i,p in zip(range(max(len(particles)-3,0),len(particles)), particles[-3:]) if int(p.type) <= 1000000001]
    # decay_products = Vector()

    lepton_length = Vector{Float64}()
    n_secondaries = Vector{Int64}()
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
            # child_particle = Particle(sec.type, sec.energy, particle, nothing)
            push!(particle.children, sec)
            # push!(decay_products,child_particle)
        end

    end
    
    E = Dict("continuous" => continuous, "epair" => epair, "brems" => brems, "ioniz" => ioniz, "photo" => photo)

    return E
end


function propagate(particle::Particle, energy::Float64, iterations::Int64, medium = nothing, propagator = nothing)

    # medium = [(1,0.5e19),(0,0.5e19),(1,0.928e19),(0,0.01289e19),(0,0.1e19)]

    # particle = Particle(1.0,1.0)

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
    # lepton.energy = 1e7
    lepton.propagated_distance = 0.0
    lepton.time = 0.0

    # make this into a dict

    lepton_length = Vector{Float64}()

    lepton.energy = energy

    # local E

    println("\n\nSTARTING PROPAGATION\n\n")

    for i in tqdm(1:1:iterations)

        # secondaries = prop.propagate(lepton).particles
        # push!(lepton_length, secondaries[end].position.magnitude()/100.0)
        # push!(n_secondaries, size(secondaries, 1))

        #lepton.energy = e
        # check if the secondaries are long-lived
        secondaries = prop.propagate(lepton)

        # logic to propagate secondaries, if necessary

        push!(lepton_length, secondaries.particles[end].position.magnitude())
        
        # global E
        # global decay_products
        E = analyze_secondaries(secondaries, particle)

        particle.E = E

    end

    event = Dict("range" => lepton_length, "E" => E, "Decay Products" => decay_products)

    return lepton_length, E, decay_products
end

end # module


# using .ProposalInterface
# using Statistics
# using DelimitedFiles

# particle = Particle(14,1.0)
# medium = [(1,1e7)]
# energy = collect(6:0.5:11)
# prop = make_propagator(particle,medium)

# lepton_length, E, decay_products = ProposalInterface.propagate(particle, 1e7, 1, medium, prop)

# print(decay_products[end].type)
























##############################################################################################################
# Muon Range Code

using .ProposalInterface
using Statistics
using DelimitedFiles

particle = Particle(14,1.0, nothing, nothing, [])

# test
dims = 3
particulas = Vector{Particle}(undef, dims)

for x in 1:dims
    particulas[x] = Particle(14,1.0, nothing, nothing, [])
end

println(particulas)

println(particulas[1] === particulas[2])
    

readline()
medium = [(0,1e7),(1,1e7)]
energy = collect(6:0.5:11)

propagate(particle, 10.0^10, 1, medium)

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