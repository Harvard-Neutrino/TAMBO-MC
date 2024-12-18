using Pkg
Pkg.activate(ENV["TAMBOSIM_PATH"] * "/Tambo")
using Tambo
using JLD2
using Glob
using ArgParse
using TOML


include("trigger_defs.jl")

interpolated_effs = load(ENV["TAMBOSIM_PATH"] * "/resources/detector_efficiencies/initial_IceTop_panel_interpolations.jld2")
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
        "--simset"
            help = "Simulation set ID"
            arg_type = String
            #required = true
        "--subsimset"
            help = "Sub-simulation set ID"
            arg_type = String
            #required = true
        "--outdir"
            help = "where to store the output"
            arg_type = String
            #default = "./"
            #required = true
        "--trigger_config"
            help = "Config file"
            arg_type = String
            default = nothing
        "--trigger_type"
            help = "Type of trigger to use"
            arg_type = String
        "--module_threshold"
            help = "Module threshold"
            arg_type = Int64
        "--event_threshold"
            help = "Event threshold"
            arg_type = Int64
    end
    return parse_args(s)
end

function load_config(file_path::String)
    if isfile(file_path)
        config_dict = TOML.parsefile(file_path)
        validate_config_file(config_dict)
        return config_dict
    else
        error("Config file not found at: $file_path")
    end
end

function setup_configuration(args)
    expected_arguments = ["simfile", "eventdictdir", "outdir", "trigger_config", "trigger_type", "module_threshold", "event_threshold"]

    config_params = Dict()

    # Check that trigger configuration is provided either via CLI or config file
    if all(isnothing, [args["trigger_config"], any(isnothing, [args["trigger_type"], args["module_threshold"], args["event_threshold"]])])
        error("Missing required arguments. Please define the trigger either via CLI or config file.")
    end

    # Check that trigger settings set by either CLI or config file, not both
    if !isnothing(args["trigger_config"]) && any(!isnothing, [args["trigger_type"], args["module_threshold"], args["event_threshold"]])
        error("Both CLI and config file interface used. Please use only one.")
    end

    # Load from config files if provided
    if !isnothing(args["trigger_config"])
        config_params = load_config(args["trigger_config"])
        for (k, v) in config_params
            if k ∉ expected_arguments
                error("Unexpected key in trigger_config file: $k")
            end
            args[k] = v
        end
    end
    return args
end

function validate_config_file(config::Dict{String, Any})
    # Check that only expected configuration parameters are present
    # so user doesn't think they're setting parameters they aren't
    expected_keys = Set(["trigger_type", "module_threshold", "event_threshold"])
    unexpected_keys = setdiff(Set(keys(config)), expected_keys)
    if !isempty(unexpected_keys)
        error("Unexpected keys found in config file: ", unexpected_keys)
    end
end

function main()
    # Parse config file and command line arguments
    args = parse_commandline()
    args = setup_configuration(args)

    sim_file = args["simfile"]
    event_dicts_path = args["eventdictdir"]
    simset_ID = args["simset"]
    subsimset_ID = args["subsimset"]
    outdir = args["outdir"]
    trigger_type = args["trigger_type"]
    module_thresh = args["module_threshold"]
    event_thresh = args["event_threshold"]

    println("Configuration:")
    println("sim_file: $sim_file")
    println("event_dicts_path: $event_dicts_path")
    println("outdir: $outdir")
    println("trigger_type: $trigger_type")
    println()

    
    event_dicts = Dict{String, Any}()
    try # TODO: bit of a hack. Should just require user to provide full event dict
        event_dicts = load(event_dicts_path * "/event_dicts_$(simset_ID)_$(subsimset_ID)_full.jld2")["hit_map"]
    catch ArgumentError
        println("No full event dict found. Falling back to old behavior of combining all event dicts here.")
        
        include("../1_make_event_dicts/combine_event_dict_files.jl")
        event_dict_files = get_filename_list(event_dicts_path, simset_ID, subsimset_ID)
        for event_dict_file in event_dict_files
            merge!(event_dicts, load(event_dict_file)["hit_map"])
        end
    end

    triggered_event_ids = []
    for (key, value) in event_dicts
        if did_trigger(value, module_thresh, event_thresh, trigger_type)
            push!(triggered_event_ids, parse(Int, split(key, "/")[end]))
        end
    end

    triggered_events = load(sim_file)["injected_events"][triggered_event_ids]

    jldopen("$(outdir)/triggered_events_$(simset_ID)_$(subsimset_ID).jld2", "w") do f
        f["triggered_events"] = triggered_events
        f["trigger_config"] = Dict(
            "module_thresh" => module_thresh,
            "event_thresh" => event_thresh,
            "trigger_type" => trigger_type
        )
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end