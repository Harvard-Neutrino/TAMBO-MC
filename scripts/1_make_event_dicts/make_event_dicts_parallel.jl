using Pkg
Pkg.activate("/n/home02/thomwg11/tambo/TAMBO-MC/Tambo")
using Tambo
using PyCall
using JLD2
using StaticArrays
using ArgParse
using Glob
using Distributions
using LinearAlgebra
using Parquet2
using DataFrames 
using ProgressBars
using TOML
using Rotations

#const zcorsika = whitepaper_normal_vec.proj
#as defined in corsika 
#const xcorsika = SVector{3}([0,-zcorsika[3]/sqrt(zcorsika[2]^2+zcorsika[3]^2),zcorsika[2]/sqrt(zcorsika[2]^2+zcorsika[3]^2)])
#const ycorsika = cross(zcorsika,xcorsika)

#const xyzcorsika = inv([
#    xcorsika.x xcorsika.y xcorsika.z;
#    ycorsika.x ycorsika.y ycorsika.z;
#    zcorsika.x zcorsika.y zcorsika.z;
#])

#const sim = jldopen("/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.jld2")
#const config = SimulationConfig(; Dict(k=>v for (k, v) in sim["config"] if k != :geo_spline_path)...)
#const geo = Tambo.Geometry("/Users/pavelzhelnin/Documents/physics/TAMBO/resources/tambo_spline.jld2", Tambo.whitepaper_coord)
#const plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)
#const altmin = 1.8925255158436627units.km
#const altmax = 4.092525515843662units.km

include("utils.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--basedir"
            help = "Directory with all the simulation directories inside"
            arg_type = String
            required = true
            #default = "/Users/pavelzhelnin/Documents/physics/TAMBO/resources/airshowers/GraphNet_00000_00001"
        "--configfile"
            help = "Injection config file"
            arg_type = String
            required = true
        "--simfile"
            help = "Simulation file, if no simulation file produced you need to specify plane coord, vector and spline path"
            arg_type = String
            required = true 
            #default = "/Users/pavelzhelnin/Documents/physics/TAMBO/resources/example_sim_file.jld2"
        "--outfile"
            help = "where to store the output"
            arg_type = String
            required = true
        "--deltas"
            help = "Distance between modules on hexagonal array"
            arg_type = Float64
            required = true
        "--length"
            help = "Length of the full array"
            arg_type = Float64
            required = true
        "--nparallel"
            help = "Number of parallel jobs happening"
            arg_type = Int
            required = true
            #default=1
        "--njob"
            help = "Which parallel job"
            arg_type = Int
            required = true
            #default=1
        "--run_desc"
            help = "A description that will become the group name for the file"
            arg_type = String
            default = ""
        "--altmin"
            help = "Minimum altitude in m"
            arg_type = Float64
            required = true
            #default = 1892.5255158436627
        "--altmax"
            help = "Maximum altitude in m"
            arg_type = Float64
            required = true
            #default = 4092.525515843662 
        "--simset"
            help = "Simset"
            arg_type = String
            required = true
        "--subsimset"
            help = "Subsimset"
            arg_type = String
            required = true
        "--size"
            help = "size of modules"
            arg_type = String 
            default = "normal"
        "--config"
            help = "Path to a TOML config file"
            arg_type = String
            default = nothing
            default = 2300
       
    end
    return parse_args(s)
end

function load_config(file_path::String)
    if isfile(file_path)
        return TOML.parsefile(file_path)
    else
        error("Config file not found at: $file_path")
    end
end

function add_hits!(d::Dict, og_df::DataFrame, modules)
    rmax = maximum([norm(m.extent) for m in modules])
    for m in modules
        df = copy(og_df)
       
        df.x .-= m.pos[1]
        df.y .-= m.pos[2]
        df.z .-= m.pos[3]

        filter!([:x,:y,:z] => (x,y,z) -> abs(x) <= rmax && abs(y) <= rmax && abs(z) <=rmax, df)
        
        if isempty(df)
            continue 
        end

        for particle in eachrow(df)
            particle[[:x,:y,:z]] = m.rot * Vector(particle[[:x,:y,:z]])
        end 

        filter!([:x,:y,:z] => (x,y,z) -> all(abs.([x,y,z]) .< m.extent / 2), df)

        if isempty(df)
            continue 
        else 
            if !(m.idx in keys(d))
                d[m.idx] = Tambo.CorsikaEvent[]
            end
            for particle in eachrow(df)
                position = SVector{3}([particle.x,particle.y,particle.z])
                push!(d[m.idx], Tambo.CorsikaEvent(particle.pdg,
                particle.kinetic_energy,position,particle.time,particle.weight))
            end 
        end 
    end 
end 


function make_hit_map(
    simset,
    subsimset,
    event_number,
    modules,
    basedir,
    xyzcorsika;
    xmax=nothing,
    xmin=nothing,
    ymin=nothing,
    ymax=nothing
)
   
    if isnothing(xmin)
        xmin = minimum([m.pos.x for m in modules])
    end
    if isnothing(xmax)
        xmax = maximum([m.pos.x for m in modules])
    end
    if isnothing(ymin)
        ymin = minimum([m.pos.y for m in modules])
    end
    if isnothing(ymax)
        ymax = maximum([m.pos.y for m in modules])
    end

    d = Dict{Int, Vector{Tambo.CorsikaEvent}}()
    files = find_extant_files(simset, subsimset, event_number, basedir)
    #for file in tqdm(files)
    for file in files
        #CORSIKA8 sometimes doesn't finish b/c job times out 
        df = DataFrame()
        try
            df = DataFrame(Parquet2.Dataset(file))
        catch
            println("This file can't be read (maybe job timed out): $file")
            continue
        end

        df = loadcorsika(select(df,Not("shower","nx","ny","nz")),xyzcorsika)
        df = filter(
                e -> xmin < e.x && e.x < xmax && ymin < e.y && e.y < ymax,
                df
            )
        add_hits!(d, df, modules)
      

        #if necessary to avoid runaway RAM usage
        #GC.gc()
        #close(pqf)
    end 
    
    return d
end

function main()
    # Parse config file and command line arguments
    # Overwrite config file values with command line arguments if provided
    expected_arguments = ["basedir", "configfile", "simfile", "outfile", "deltas", "length", "nparallel", "njob", "run_desc", "altmin", "altmax", "simset", "subsimset", "config"]
    args = parse_commandline()
    
    config_params = Dict()

    # Load from config files if provided
    if haskey(args, "config") && !isnothing(args["config"])
        println(args)
        config_params = load_config(args["config"])
        for (k, v) in config_params
            if k ∉ expected_arguments
                error("Unexpected key in config file: $k")
            end
            args[k] = v
        end
    end

    # Override with command line arguments if provided
    for (k, v) in args
        if k ∉ expected_arguments
            error("Unexpected key in command line arguments: $k")
        end
        if k != "config"
            config_params[k] = v
        end
    end
    
    @assert args["njob"] <= args["nparallel"]

    println("Configuration:")
    for (k, v) in args
        println("$k: $v")
    end

    #should fix so that we have a .toml 
    sim = jldopen(args["simfile"])
    config = Simulation("/n/home02/thomwg11/tambo/TAMBO-MC/resources/configuration_examples/$(args["configfile"]).toml")
    geo = Tambo.Geometry(config.config["geometry"])
    tambo_coord_degrees = Tambo.Coord((deg2rad.(config.config["geometry"]["tambo_coordinates"]))...)
    plane = Tambo.Plane(Tambo.Direction(config.config["geometry"]["plane_orientation"]...), tambo_coord_degrees, geo)
    altmin = args["altmin"]units.m
    altmax = args["altmax"]units.m
    
    zcorsika = config.config["geometry"]["plane_orientation"]
    display("original zcorsika = $zcorsika")
    #possibly a hack but (this is because everything has to be rotated to corsika x_hat = N;y_hat=W)
    #previous convention was x_hat = E; y_hat =N 
    #xyzcorsika rotates into 3d from 2d but to where x_hat is N!!!!
    #need extra rotation to get it back to TAMBO coordinate system 
    
    zcorsika = RotZ(-π/2) * zcorsika 
    xcorsika = SVector{3}([0,-zcorsika[3]/sqrt(zcorsika[2]^2+zcorsika[3]^2),zcorsika[2]/sqrt(zcorsika[2]^2+zcorsika[3]^2)])
    ycorsika = cross(zcorsika,xcorsika)
    xyzcorsika = inv([
        xcorsika.x xcorsika.y xcorsika.z;
        ycorsika.x ycorsika.y ycorsika.z;
        zcorsika.x zcorsika.y zcorsika.z;
    ])

    display("xyzcorsika = $xyzcorsika")

    if args["size"] == "normal"
        size = SVector{3}([1.875, 0.8, 0.03])units.m
    elseif args["size"] == "small"
        size = SVector{3}([2,2, 0.03])units.m
    elseif args["size"] == "medium"
        size = SVector{3}([4,4,0.03])units.m
    end 

    modules = Tambo.make_detector_array(
        args["length"]units.m,
        args["deltas"]units.m,
        altmin,
        altmax,
        plane,
        geo,
        size,
    )

    outfile = args["outfile"]
    #if args["nparallel"] > 1
    #    outfile = replace(outfile, ".jld2"=>"_$(args["njob"])_$(args["nparallel"]).jld2")
    #end
    
    run_desc = args["run_desc"]
    if length(run_desc)==0
        run_desc = "$(args["length"])_$(args["deltas"])"
    end

    event_numbers = get_event_numbers(args["basedir"], args["simset"], args["subsimset"])[args["njob"]:args["nparallel"]:end]
    println("creating event dicts...")
    println("event numbers: $event_numbers")

    hit_map = Dict()
    #for (_,event_number) in tqdm(enumerate(event_numbers))
    for event_number in event_numbers
        hit_map["$(run_desc)/$(event_number)"] = make_hit_map(args["simset"], args["subsimset"], event_number, modules, args["basedir"], xyzcorsika)
    end
    jldopen(outfile, "w") do jldf
        for (k, v) in hit_map
            jldf[k] = v
        end

    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
