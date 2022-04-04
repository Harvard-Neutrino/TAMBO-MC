using PyCall
using ProgressBars
using DelimitedFiles

pp = pyimport("proposal")


# Make this useful
struct Particle
    pdg_mc::Float64
    energy::Float64
end

function make_sector(density::Float64, start::Float64, stop::Float64)

    # Defining a Sector
    sec_def = pp.SectorDefinition()

    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), stop, start)

    # I have no idea what this does
    components = [pp.component.Hydrogen(2),pp.component.Oxygen()]
    sec_def.medium = pp.medium.Medium("Ice_$density", 1, 75.0, -3.5017, 
    0.09116, 3.4773, 0.2400, 
    2.8004,0, density, components)

    # sec_def.medium = pp.medium.Ice(1.0)

    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere

    # The way it was
    sec_def.crosssection_defs.brems_def.lpm_effect = true
    sec_def.crosssection_defs.epair_def.lpm_effect = true
    sec_def.cut_settings.ecut = 0
    sec_def.cut_settings.vcut = 0.05
    # sec_def.crosssection_defs.brems_def.lpm_effect = false
    # sec_def.crosssection_defs.epair_def.lpm_effect = false
    # sec_def.cut_settings.ecut = 500
    # sec_def.cut_settings.vcut = 0.05

    sec_def.do_continuous_randomization = true

    return sec_def

end

# Make this useful
function define_particle(particle::Particle)
    return pp.particle.MuMinusDef()
end


function make_propagator(particle::Particle, medium::Vector{Tuple{Float64,Float64}})

    # Definig sectors from medium (Array) and calculating detector length
    sectors = Vector{PyObject}()

    println("Defining Sectors:")

    detector_length = 0.0
    sector_count = 0
    for (ρ,l) in medium

        sector_count += 1
        # detector_length
        println("Sector $sector_count")
        println("Density: $ρ")
        println("Length: $l")

        # Calculating start and stop points for each sector
        start = detector_length
        # updates detector length
        detector_length += l
        stop = detector_length
        push!(sectors,make_sector(ρ,start,stop))
    end

    #println(sectors)
    println("Detector Length: $detector_length")

    particle = Particle(1.0,1.0)

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


# Testing it out plotting MuMinus ranges

medium = [(0.97,0.5e19),(1.19,0.5e19)]

particle = Particle(1.0,1.0)

particle_def = define_particle(particle)

prop = make_propagator(particle,medium)


mu = pp.particle.DynamicData(particle_def.particle_type)
mu.position = pp.Vector3D(0, 0, 0)
mu.direction = pp.Vector3D(0, 0, -1)
mu.energy = 1e7
mu.propagated_distance = 0.0
mu.time = 0.0

# mu_energies = Vector{Float64}()
mu_length = Vector{Float64}()
n_secondaries = Vector{Int64}()

# pp.RandomGenerator.get().set_seed(1234)

println("\n\nSTARTING PROPAGATION\n\n")

for i in tqdm(1:10000)

    secondaries = prop.propagate(mu).particles
    push!(mu_length, secondaries[end].position.magnitude()/100.0)
    push!(n_secondaries, size(secondaries, 1))

    # secondaries = prop.propagate(mu,30000)
    # secondaries = prop.propagate(mu,3000)
    # energy = last(secondaries.energy)
    # # println(energy)
    # push!(mu_energies, energy)

end

# print(mu_length)
# print(n_secondaries)

writedlm( "Results/MuMinus_Length.csv",  mu_length, ',')
writedlm( "Results/MuMinus_Secondaries.csv",  n_secondaries, ',')

println("Files saved in Results Folder")


