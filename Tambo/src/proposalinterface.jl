using PyCall
#using ProgressBars
using DelimitedFiles

#export Particle, make_sector, define_particle, make_medium, make_propagator, propagate

pp = pyimport("proposal")


const propagators = Dict()

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
    sec_def.crosssection_defs.photo_def.parametrization = (
        pp.parametrization.photonuclear.PhotoParametrization.AbramowiczLevinLevyMaor97
    )
    sec_def
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
    eval(particles[particle.pdg_mc])
end


function make_medium(medium)

    # Defining sectors from medium (Array) and calculating detector length
    sectors = Vector()

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

    losses = [[], [], [], [], []]

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
                child = Particle(sec.type + 1, sec.energy, 0.0, 0.0, parent_particle, [])
                push!(parent_particle.children, child)
            elseif sec.type in [-13, -15, -11] 
                child = Particle(sec.type - 1, sec.energy, 0.0, 0.0, parent_particle, [])
                push!(parent_particle.children, child)
            end
        else
            println("$(sec.type)")
        end
        

    end

    parent_particle.final_vertex = SVector{3}(secondaries.particles[end].position)

    E = Dict("continuous" => continuous, "epair" => epair, "brems" => brems, "ioniz" => ioniz, "photo" => photo)

    parent_particle.E = E

end


function propagate(particle::Particle, medium = nothing, propagator = nothing)

    particle_def = define_particle(particle)

    # Stores previous propagators to make code faster
    if particle.pdg_mc in keys(propagators)
        prop = propagators[particle.pdg_mc]
    else
        prop = make_propagator(particle,medium)
        propagators[particle.pdg_mc] = prop
    end

    lepton = pp.particle.DynamicData(particle_def.particle_type)
    lepton.position = pp.Vector3D(0, 0, 0)
    lepton.direction = pp.Vector3D(0, 0, -1)
    lepton.energy = particle.energy
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
        "Initial Energy" : $(particle.energy),
        "Range" : $(particle.range),
        "Energy breakdown" : $(particle.E),
        "Decay Products" : $(particle.children))
        }""")
end