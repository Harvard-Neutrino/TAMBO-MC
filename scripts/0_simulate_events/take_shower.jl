using Pkg
Pkg.activate(".")
Pkg.develop(path=ENV["TAMBOSIM_PATH"] * "/Tambo")
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
        "--simset_id"
            help = "Simulation set ID"
            arg_type = Int
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
    simset_id = args["simset_id"]

    sim = Simulation(config_filename, injection_filename)

    # A single pinecone is used to generate a unique seed for each simulation
    # based on the simulation set ID and sub-simulation set ID.
    # It is a pinecone because pinecones release seeds.
    pinecone = sim.config["steering"]["pinecone"]

    seed!(pinecone)
    seed = rand(1:typemax(Int64)) + simset_id
    seed!(seed)

    sim.config["corsika"]["shower_dir"] = shower_dir

    run_subshower!(sim, sim.config["corsika"], proposal_id, decay_id, seed)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
