using Pkg
Pkg.activate("../../Tambo")
using Tambo
using Glob
using JLD2
include("trigger_defs.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--basedir"
            help = "Directory with all the simulation directories inside"
            arg_type = String
            #required = true
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/event_dicts"
        "--eventdictdir"
            help = "Directory with all the event dicts inside"
            arg_type = String
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/sim_files"
            #required = true 
        "--outdir"
            help = "where to store the output"
            arg_type = String
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/triggered_events"
            #required = true
    end
    return parse_args(s)
end

function main()
    outdir = args["outfile"]

    if !isdir(outdir)
        # Create the directory if it does not exist
        mkpath(outdir)
    end

    files = glob("*.jld2",args["eventdictdir"])
    
    triggered_events = Dict{String, Vector{Int}}()
    for file in files 
        file_number = split(file,"_")[end]
        f = jldopen(file)
        for det_config in keys(f)
            for event in keys(f[det_config])
                if did_trigger(f[det_config][event])
                    if !(file_number in keys(triggered_events))
                        triggered_events[file_number] = []
                    end
                    push!(triggered_events[file_number],Int64(event))
                end
            end
        end
    end 
    println(triggered_events)
end
