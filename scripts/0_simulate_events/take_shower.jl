using Pkg
Pkg.activate("/n/home02/thomwg11/tambo/TAMBO-MC/Tambo")
using Tambo
using ArgParse
using Random: seed!

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config"
            help = "Path to a TOML config file"
            arg_type = String
            required = true
        "--injection_file"
            help = "Path to a JLD2 file containing the injection data"
            arg_type = String
            required = true
        "--shower_dir"
            help = "Path to directory in which to save the output"
            arg_type = String
            required = true
        "--proposal_id"
            help = "ID of the tau from PROPOSAL to simulate"
            arg_type = Int
            required = true
        "--decay_id"
            help = "ID of the decay particle from the tau to simulate"
            arg_type = Int
            required = true
        "--simset"
            help = "Simulation set ID"
            arg_type = String
            required = true
        "--subsimset"
            help = "Sub-simulation set ID"
            arg_type = String
            required = true
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
    args = parse_commandline()
    config_filename = args["config"]
    injection_filename = args["injection_file"]
    shower_dir = args["shower_dir"]
    proposal_id = args["proposal_id"]
    decay_id = args["decay_id"]
    simset_ID = args["simset"]
    subsimset_ID = args["subsimset"]

    sim = Simulation(config_filename, injection_filename)

    seed = sim.config["steering"]["seed"]
    seed!(seed + parse(Int,simset_ID) + parse(Int,subsimset_ID)) # TODO: not a unique seed

    # Hack to make sure user doesn't try to set output_dir in config file
    # TODO: remove support for setting output_dir in config file
    if sim.config["corsika"]["shower_dir"] != ""::String
        error("output_dir must be set via command line argument, not in the config file")
    end

    sim.config["corsika"]["shower_dir"] = shower_dir

    run_subshower!(sim, sim.config["corsika"], proposal_id, decay_id)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end