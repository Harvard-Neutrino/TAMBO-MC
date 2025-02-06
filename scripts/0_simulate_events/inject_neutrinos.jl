using Pkg
Pkg.activate(ENV["TAMBOSIM_PATH"] * "/scripts/0_simulate_events/")
Pkg.develop(path=ENV["TAMBOSIM_PATH"] * "/Tambo")
Pkg.add("ArgParse")
#Pkg.activate(ENV["TAMBOSIM_PATH"] * "/Tambo")
using Tambo
using ArgParse
using Random: seed!

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config"
            help = "Path to a TOML config file"
            arg_type = String
            default = nothing
        "--simset"
            help = "Simulation set ID"
            arg_type = Int
            required = true
        "--subsimset"
            help = "Sub-simulation set ID"
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
    #if match(r"^(.+)_\d{5}_\d{5}\((.jld2)|(.arrow))", output_filename) === nothing
    if match(r"^(.+)_\d{5}_\d{5}\.((jld2)|(arrow))$", output_filename) === nothing
        error("Output filename must adhere to format <name>_xxxxx_yyyyy.jld2")
    end
end

function main()
    args = parse_commandline()
    config_filename = args["config"]
    simset_ID = args["simset"]
    subsimset_ID = args["subsimset"]
    output_filename = args["output"]

    validate_output_filename(output_filename)

    sim = Simulation(config_filename)

    # A single pinecone is used to generate a unique seed for each simulation
    # based on the simulation set ID and sub-simulation set ID.
    # It is a pinecone because pinecones release seeds.
    pinecone = sim.config["steering"]["pinecone"]

    seed!(pinecone)
    seed = round(Int64, rand()) + 100000*simset_ID + subsimset_ID
    seed!(seed)

    inject_ν!(sim, sim.config["injection"], seed; track_progress=true)
    propagate_τ!(sim, sim.config["proposal"], seed)
    identify_taus_to_shower!(sim, sim.config["corsika"], seed)

    # Create output directory if it does not exist
    output_dir = dirname(output_filename)
    if output_dir != "" && ~isdir(output_dir)
        mkdir(output_dir)
    end

    if endswith(output_filename, "jld2")
        save_simulation_to_jld2(sim, output_filename)
    elseif endswith(output_filename, "arrow")
        save_simulation_to_arrow(sim, output_filename)
    else
        save_simulation_to_jld2(sim, output_filename)
        save_simulation_to_arrow(sim, output_filename)
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
