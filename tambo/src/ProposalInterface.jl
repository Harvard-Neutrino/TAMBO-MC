module ProposalInterface

using PyCall
using ProgressBars
using DelimitedFiles

export Particle, make_sector, define_particle, make_medium, make_propagator, propagate_mcp

pp = pyimport("proposal")


# Make this useful
struct Particle
    pdg_mc::Float64
    energy::Float64
end

function make_sector(medium, start, stop)

    # Defining a Sector
    sec_def = pp.SectorDefinition()

    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), stop, start)

    if medium == 1
        # sec_def.medium = pp.medium.StandardRock()

        # Muon Test only
        sec_def.medium = pp.medium.Ice()
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

# Make this useful
function define_particle(particle::Particle)
    # {16: 'TauMinusDef', -16: 'TauPlusDef',
    #          14: 'MuMinusDef',  -14: 'MuPlusDef',
    #          12: 'EMinusDef',   -12: 'EPlusDef',}
    return pp.particle.MuMinusDef()
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
        push!(sectors,make_sector(rock,start,stop))
    end

    #println(sectors)
    println("Detector Length: $detector_length")

    return sectors, detector_length
end


function make_propagator(particle::Particle, medium)

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


function propagate_mcp(particle::Particle, energy::Float64, iterations::Int64, medium = nothing, propagator = nothing)

    # medium = [(1,0.5e19),(0,0.5e19),(1,0.928e19),(0,0.01289e19),(0,0.1e19)]

    # particle = Particle(1.0,1.0)

    particle_def = define_particle(particle)
    if propagator == nothing
        prop = make_propagator(particle,medium)
    else
        prop = propagator
    end

    # Probably gonna have to fix this.
    mu = pp.particle.DynamicData(particle_def.particle_type)
    mu.position = pp.Vector3D(0, 0, 0)
    mu.direction = pp.Vector3D(0, 0, -1)
    # mu.energy = 1e7
    mu.propagated_distance = 0.0
    mu.time = 0.0

    # mu_energies = Vector{Float64}()
    mu_length = Vector{Float64}()
    n_secondaries = Vector{Int64}()

    continuous = Vector{Vector}()
    epair = Vector{Vector}()
    brems = Vector{Vector}()
    ioniz = Vector{Vector}()
    photo = Vector{Vector}()

    # pp.RandomGenerator.get().set_seed(1234)
    mu.energy = energy

    println("\n\nSTARTING PROPAGATION\n\n")

    for i in tqdm(1:1:iterations)

        # secondaries = prop.propagate(mu).particles
        # push!(mu_length, secondaries[end].position.magnitude()/100.0)
        # push!(n_secondaries, size(secondaries, 1))

        #mu.energy = e

        secondaries = prop.propagate(mu)

        push!(mu_length, secondaries.particles[end].position.magnitude()/100.0)
        push!(n_secondaries, size(secondaries.particles, 1))

        for sec in secondaries.particles

            log_sec_energy = log.(10,sec.parent_particle_energy .- sec.energy)
            
            if sec.type == convert(Int64, pp.particle.Interaction_Type.ContinuousEnergyLoss)
                push!(continuous, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
            end
            if sec.type == convert(Int64, pp.particle.Interaction_Type.Epair)
                push!(epair, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
            end
            if sec.type == convert(Int64, pp.particle.Interaction_Type.Brems)
                push!(brems, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
            end
            if sec.type == convert(Int64, pp.particle.Interaction_Type.DeltaE)
                push!(ioniz, [log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
            end
            if sec.type == convert(Int64, pp.particle.Interaction_Type.NuclInt)
                push!(photo,[log_sec_energy, [sec.position.x, sec.position.y, sec.position.z]])
            end

        end

        # secondaries = prop.propagate(mu,30000)
        # secondaries = prop.propagate(mu,3000)
        # energy = last(secondaries.energy)
        # # println(energy)
        # push!(mu_energies, energy)

    end

    E = Dict("continuous" => continuous, "epair" => epair, "brems" => brems, "ioniz" => ioniz, "photo" => photo)

    # writedlm( "Results/MuMinus_Length.csv",  mu_length, ',')
    # writedlm( "Results/MuMinus_Secondaries.csv",  n_secondaries, ',')

    # println("Files saved in Results Folder")

    return mu_length, n_secondaries, E
end

end # module

using .ProposalInterface
using Statistics

particle = Particle(1.0,1.0)
medium = [(1,1e6)]
energy = collect(6:0.33:11)
prop = make_propagator(particle,medium)

energies = Vector{Float64}()
range = Vector{Float64}()
error = Vector{Float64}()

for e in energy
    println("\n\n\nLogEnergy: $e MeV")
    mu_length, n_secondaries, E = ProposalInterface.propagate_mcp(particle, 10^e, 500, medium, prop)
    d = mean(mu_length)
    err = std(mu_length)
    println("Mean Range: $d")
    println("Standard Deviation: $err")
    push!(range, d)
    push!(energies, e)
    push!(error, err)
end

println(error)
println(range)
println(energies)

# writedlm( "Results/MuMinus_Length.csv",  mu_length, ',')
# writedlm( "Results/MuMinus_Secondaries.csv",  n_secondaries, ',')

# println("Files saved in Results Folder")

# mu_length, n_secondaries, E = propagate_mcp(particle, medium, mcp_energies)

# print(mu_length)
# print(n_secondaries)
# print(E)