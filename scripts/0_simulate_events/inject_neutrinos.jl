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
    if match(r"^(.+)_\d{5}_\d{5}\.jld2$", output_filename) === nothing
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

    seed = sim.config["steering"]["seed"]

    seed!(seed + simset_ID + subsimset_ID) # TODO: not a unique seed

    inject_ν!(sim, sim.config["injection"])
    propagate_τ!(sim, sim.config["proposal"])
    #run_airshower!(sim, sim.config["corsika"])
    identify_taus_to_shower!(sim, sim.config["corsika"])

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
