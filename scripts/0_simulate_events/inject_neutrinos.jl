using Pkg
using Random: seed!
Pkg.activate(ENV["TAMBOSIM_PATH"] * "/Tambo")
using Tambo
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config"
            help = "Path to a TOML config file"
            arg_type = String
            default = nothing
        "--simset_id"
            help = "Simulation set ID"
            arg_type = Int
            required = true
        "--output"
            help = "Output file"
            arg_type = String
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

function validate_output_filename(output_filename::String)
    # Check that output filename adheres to specified <name>_xxxxx_yyyyy.jld2 format
    if match(r"^(.+)_\d{5}_\d{5}\.jld2$", output_filename) === nothing
        error("Output filename must adhere to format <name>_xxxxx_yyyyy.jld2")
    end
end

function main()
    args = parse_commandline()
    config_filename = args["config"]
    simset_id = args["simset_id"]
    output_filename = args["output"]

    validate_output_filename(output_filename)

    sim = Simulation(config_filename)

    # A single pinecone is used to generate a unique seed for each simulation
    # based on the simulation set ID and sub-simulation set ID.
    # It is a pinecone because pinecones release seeds.
    pinecone = sim.config["steering"]["pinecone"]

    seed = pinecone + simset_id
    seed!(seed)

    inject_ν!(sim, sim.config["injection"], seed)
    propagate_τ!(sim, sim.config["proposal"], seed)
    identify_taus_to_shower!(sim, sim.config["corsika"], seed)

    # Create output directory if it does not exist
    output_dir = dirname(output_filename)
    if output_dir != "" && ~isdir(output_dir)
        mkdir(output_dir)
    end

    save_simulation(sim, output_filename)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
