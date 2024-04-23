using Pkg
Pkg.activate("../Tambo")
using Tambo
using PyCall
using JLD2
using StaticArrays
using ArgParse
using Glob
using Distributions
using LinearAlgebra

const pq = PyNULL()
const np = PyNULL()
copy!(pq, pyimport("pyarrow.parquet"))
copy!(np, pyimport("numpy"))

const ycorsika = SVector{3}([0.89192975455881607, 0.18563051261662877, -0.41231374670066206])
const xcorsika = SVector{3}([0, -0.91184756344828699, -0.41052895273466672])
const zcorsika = whitepaper_normal_vec.proj
const xyzcorsika = inv([
    xcorsika.x xcorsika.y xcorsika.z;
    ycorsika.x ycorsika.y ycorsika.z;
    zcorsika.x zcorsika.y zcorsika.z;
])
const sim = jldopen("/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.jld2")
const config = SimulationConfig(; Dict(k=>v for (k, v) in sim["config"] if k != :geo_spline_path)...)
const geo = Tambo.Geometry(config)
const plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)

const altmin = 1.8925255158436627units.km
const altmax = 4.092525515843662units.km

include("utils.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--basedir"
            help = "Directory with all the simulation directories inside"
            arg_type = String
            required = true
        "--outfile"
            help = "where to store the output"
            arg_type = String
            required = true
        "--deltas"
            help = "Distance between modules on hexagonal array"
            arg_type = Float64
            required = true
        "--length"
            help = "Length of the full array"
            arg_type = Float64
            required = true
        "--nparallel"
            help = "Number of parallel jobs happening"
            arg_type = Int
            default=1
        "--njob"
            help = "Which parallel job"
            arg_type = Int
            default=1
        "--run_desc"
            help = "A description that will become the group name for the file"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function add_hits!(d::Dict, events, modules)
    rmax = maximum([norm(m.extent) for m in modules])
    for (idx, event) in enumerate(events)
        a = inside.(Ref(event.pos), modules, rmax)
        s = sum(a)
        @assert s <= 1 "Seems like you're in more than one module.. That doesn't seem right"
        if s > 0
            mod = modules[findfirst(a)]
            if !(mod.idx in keys(d))
                d[mod.idx] = Tambo.CorsikaEvent[]
            end
            push!(d[mod.idx], event)
        end
    end
end

function make_hit_map(
    event_number,
    modules,
    basedir;
    xmax=nothing,
    xmin=nothing,
    ymin=nothing,
    ymax=nothing
)
    if isnothing(xmin)
        xmin = minimum([m.pos.x for m in modules])
    end
    if isnothing(xmax)
        xmax = maximum([m.pos.x for m in modules])
    end
    if isnothing(ymin)
        ymin = minimum([m.pos.y for m in modules])
    end
    if isnothing(ymax)
        ymax = maximum([m.pos.y for m in modules])
    end

    d = Dict{Int, Vector{Tambo.CorsikaEvent}}()
    files = find_extant_files(event_number, basedir)
    for file in files
        # This nastiness needs to be here cuz sometimes the file is empty(?)
        # and it will error in an annoying way. I think the `pqf = nothing` is superfluous
        # but I put it there before and I'm worried it will cause issues
        pqf = nothing
        #try
            pqf = pq.ParquetFile(file)
        #catch
        #    continue
        #end
        jdx = 1
        for batch in pqf.iter_batches()
            events = loadcorsika(batch)
            jdx += 1
            events = filter(
                e -> xmin < e.pos.x && e.pos.x < xmax && ymin < e.pos.y && e.pos.y < ymax,
                events
            )
            add_hits!(d, events, modules)
            # These two lines avoid runaway RAM usage
            events = nothing
            GC.gc()
        end
    end
    return d
end

function can_skip_event(event_number, outfile, run_desc)
    can_skip = false
    jldopen(outfile, "r") do jldf
        if ~(run_desc in keys(jldf))
            can_skip = false
        elseif string(event_number) in keys(jldf[run_desc])
            can_skip = true
        end
    end
    return can_skip
end

function main()
    args = parse_commandline()
    
    @assert args["njob"] <= args["nparallel"]
    modules = Tambo.make_detector_array(
        whitepaper_coord,
        args["length"]units.m,
        args["deltas"]units.m,
        altmin,
        altmax,
        plane,
        geo
    )
    outfile = args["outfile"]
    if args["nparallel"] > 1
        outfile = replace(outfile, ".jld2"=>"_$(args["njob"])_$(args["nparallel"]).jld2")
    end

    if ~ispath(outfile)
        jldopen(outfile, "w") do _
        end
    end
    
    run_desc = args["run_desc"]
    if length(run_desc)==0
        run_desc = "$(args["length"])_$(args["deltas"])"
    end

    event_numbers = get_event_numbers(args["basedir"])[args["njob"]:args["nparallel"]:end]
    for event_number in event_numbers
        @show event_number

        if can_skip_event(event_number, outfile, run_desc)
            continue
        end
        hit_map = make_hit_map(event_number, modules, args["basedir"])
        if length(hit_map)==0
            continue
        end
        jldopen(outfile, "r+") do jldf
            jldf["$(run_desc)/$(event_number)"] = hit_map
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
