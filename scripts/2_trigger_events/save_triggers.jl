using Pkg
Pkg.activate("../../Tambo")
using Tambo
using Glob
using JLD2
include("trigger_defs.jl")
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--eventdictdir"
            help = "Directory with all the simulation directories inside"
            arg_type = String
            #required = true
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/event_dicts"
        "--simfiles"
            help = "Directory with all the event dicts inside"
            arg_type = String
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/sim_files"
            #required = true 
        "--outdir"
            help = "where to store the output"
            arg_type = String
            #default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/triggered_events"
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/TAMBO-MC/scripts/test_triggered_events"
            #required = true
    end
    return parse_args(s)
end

function make_triggered_event_dicts(d::Dict,simfiles)
for det_config in keys(d)
    for ending in keys(d[det_config])
        events = Vector{Tambo.InjectionEvent}()
        filename = joinpath(simfiles,"larger_valley_00000_$(ending)")
        file = jldopen(filename)
        for event in d[det_config][ending]
            push!(events, file["injected_events"][event])
        end
    # Save the dictionary to a JLD2 file
        ending = split(ending,".")[1]
        jldopen("larger_valley_00000_$(ending)_triggered_events.jld2", "w") do file
            write(file, det_config, events)
        end
    end
end 

end

function main()
    args = parse_commandline()
    outdir = args["outdir"]

    if !isdir(outdir)
        # Create the directory if it does not exist
        mkdir(outdir)
    end

    files = glob("larger*10.jld2",args["eventdictdir"])
    simfiles = glob(".jld2",args["simfiles"])
    #triggered_events = Dict{String, Vector{Int}}()
    triggered_events = Dict{String, Dict{String, Vector{Int}}}()
    for file in files 
        println(file)
        file_number = split(file,"_")[end]
        f = jldopen(file)

        #multiple detector configs in a single event dict file 
        for det_config in keys(f)
            #create a dict if det config not in triggered_events
            if !(det_config in keys(triggered_events))
                triggered_events[det_config] = Dict{String,Vector{Int}}()
            end
            
            #looping through events in each det_config 
            for event in keys(f[det_config])

                if did_trigger(f[det_config][event])
                    if !(file_number in keys(triggered_events[det_config]))
                        triggered_events[det_config][file_number] = []
                    end
                    push!(triggered_events[det_config][file_number],parse(Int64,event))
                end
            end
        end
    end 
    println(triggered_events)
    make_triggered_event_dicts(triggered_events,simfiles)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
