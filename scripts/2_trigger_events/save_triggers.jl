using Pkg
Pkg.activate("../../Tambo")
using Tambo
using Glob
using JLD2
using ArgParse
using ProgressBars
include("trigger_defs.jl")
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--eventdictdir"
            help = "Directory with all the simulation directories inside"
            arg_type = String
            #required = true
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/near_hard/super_long_event_dicts"
        "--basedir"
            help = "Directory with all the event dicts inside"
            arg_type = String
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/near_hard/sim_files"
            #required = true 
        "--outdir"
            help = "where to store the output"
            arg_type = String
            #default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/triggered_events"
            default = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/near_hard/triggered_files"
            #required = true
    end
    return parse_args(s)
end

function return_events(files)
    triggered_events = Dict{String, Vector{Int}}()
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
                    push!(triggered_events[file_number],parse(Int,event))
                end
            end
        end
    end 
    return triggered_events
end 

function make_triggered_event_dicts(d::Dict,outfile::String,simfile::String)
    println(keys(d))
    for ending in keys(d["5000.0_150.0"])
        events = Vector{Tambo.InjectionEvent}()
        filename = joinpath(simfile,"larger_valley_00000_$(ending)")
        file = jldopen(filename)
        for event in d["5000.0_150.0"][ending]
            push!(events, file["injected_events"][event])
        end
    # Save the dictionary to a JLD2 file
        ending = split(ending,".")[1]
        jldopen("$(outfile)/retest_normal_larger_valley_00000_$(ending)_triggered_events.jld2", "w") do file
            write(file, "5000.0_150.0", events)
        end
    end
end

function main()
    args = parse_commandline()
    simfile = args["basedir"]
    outfile = args["outdir"]
    if !isdir(outfile)
        # Create the directory if it does not exist
        mkpath(outfile)
    end


    println(args["eventdictdir"])
    files = glob("retest_normal_larger_valley_event_dicts_00000_000*.jld2",args["eventdictdir"])
    triggered_events = Dict{String, Dict{String, Vector{Int}}}()
    #simfile = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/pavel_cross_check/"
    #files = ["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Larger_Valley/event_dicts/Will_normal_larger_valley_event_dicts_00000_00001.jld2"]
    
    for file in tqdm(files) 
        file_number = split(file,"_")[end]                                                                              
        f = jldopen(file)

        for det_config in keys(f)
            if !(det_config in keys(triggered_events))
                triggered_events[det_config] = Dict{String,Vector{Int}}()
            end 
            for event in keys(f[det_config])
                
                #println(did_trigger(f[det_config][event]))
            
                if did_trigger(f[det_config][event],3,30)
                    if !(file_number in keys(triggered_events[det_config]))
                        triggered_events[det_config][file_number] = []
                    end
                    push!(triggered_events[det_config][file_number],parse(Int64,event))
                end
            end
        end
    end 
    println(triggered_events)
    make_triggered_event_dicts(triggered_events,outfile,simfile)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
