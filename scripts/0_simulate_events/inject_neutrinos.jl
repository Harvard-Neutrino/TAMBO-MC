using Pkg
Pkg.activate(ENV["TAMBOSIM_PATH"] * "/scripts/0_simulate_events/")
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
    # Check that output filename adheres to specified <name>_xxxxx.jld2 format
    if match(r"^(.+)_\d{5}\.(jld2|arrow)$", output_filename) === nothing
        error("Output filename must adhere to format <name>_xxxxx.jld2")
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

    # Don't really want our different seeds to be right next to each otherwise
    # we might get correlated results. Space out the seeds by seeding the RNG
    # with pinecone + simset_id first and then drawing a random seed from that.
    seed!(pinecone + simset_id)
    # Max seed value is typemax(Int32) so we subtract 1_000_000 to avoid overflow.
    # I'm assuming 1_000_000 is the largest value for the simset_id we'll ever use.
    seed = rand(0:typemax(Int32)) - 1_000_000 + simset_id 

    seed!(seed)
    inject_ν!(sim, sim.config["injection"], simset_id, seed)
    seed!(seed) # reseed to prevent dependence on number of events injected
    propagate_τ!(sim, sim.config["proposal"], seed)
    seed!(seed) # reseed to prevent dependence on number of taus propagated
    identify_taus_to_shower!(sim, sim.config["corsika"])

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
