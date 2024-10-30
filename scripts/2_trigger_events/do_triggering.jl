using Pkg
Pkg.activate("../Tambo")
using Tambo
using JLD2
using Glob
using ArgParse
using TOML


include("trigger_defs.jl")

interpolated_effs = load("../resources/detector_efficiencies/initial_IceTop_panel_interpolations.jld2")
global interpolated_eff_gamma = interpolated_effs["gamma_interp"]
global interpolated_eff_muon = interpolated_effs["muon_interp"]
global interpolated_eff_electron = interpolated_effs["electron_interp"]

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--simfile"
            help = "Main JLD2 simulation file"
            arg_type = String
            #required = true
            #default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Oct16th2023_WhitePaper_300k.jld2"
        "--eventdictdir"
            help = "Directory with all the event dicts inside"
            arg_type = String
            #default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/triggered_events/Jan7th2024_WhitePaper_300k_no_thin/"
            #required = true 
        "--outdir"
            help = "where to store the output"
            arg_type = String
            #default = "./"
            #required = true
        "--trigger_type"
            help = "Type of trigger to use"
            arg_type = String
            #default = "whitepaper"
        "--config"
            help = "Config file"
            arg_type = String
            default = nothing
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

function main()
    # Parse config file and command line arguments
    # Overwrite config file values with command line arguments if provided
    expected_arguments = ["simfile", "eventdictdir", "outdir", "trigger_type", "config"]
    args = parse_commandline()
    
    config_params = Dict()

    # Load from config files if provided
    if haskey(args, "config") && !isnothing(args["config"])
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

    sim_file = args["simfile"]
    event_dicts_path = args["eventdictdir"]
    outdir = args["outdir"]
    trigger_type = args["trigger_type"]

    println("Configuration:")
    println("sim_file: $sim_file")
    println("event_dicts_path: $event_dicts_path")
    println("outdir: $outdir")
    println("trigger_type: $trigger_type")
    println()

    if trigger_type == "whitepaper"
        module_trigger_thresh = 3
        event_trigger_thresh = 30
    
    elseif trigger_type == "threeamigos"
        module_trigger_thresh = 3
        event_trigger_thresh = 30

    elseif trigger_type == "icetop_tanks"
        module_trigger_thresh = 300
        event_trigger_thresh = 3000

    elseif trigger_type == "icetop_panels"
        module_trigger_thresh = 3 * 65
        event_trigger_thresh = 30 * 65
    end

    sim = jldopen(sim_file)
    
    event_dict_files = glob("*.jld2", event_dicts_path)
    event_dicts = Dict{String, Any}()


    for event_dict_file in event_dict_files
        merge!(event_dicts, load(event_dict_file))
    end

    triggered_event_ids = []
    for (key, value) in event_dicts
        if did_trigger(value, module_trigger_thresh, event_trigger_thresh, trigger_type)
            push!(triggered_event_ids, parse(Int, split(key, "/")[end]))
        end
    end
    triggered_events = sim["injected_events"][triggered_event_ids]

    save("$(outdir)/triggered_events_$(trigger_type)_mod_$(module_trigger_thresh)_event_$(event_trigger_thresh).jld2", "triggered_events", triggered_events)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end