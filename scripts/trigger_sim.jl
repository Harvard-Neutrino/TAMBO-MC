using Pkg
Pkg.activate("../Tambo")
using Tambo
using PyCall
using JLD2
using StaticArrays
using ArgParse
using Glob
using Distributions

const pq = PyNULL()
const np = PyNULL()
copy!(pq, pyimport("pyarrow.parquet"))
copy!(np, pyimport("numpy"))

const zmin = -1100units.m
const zmax = 1100units.m
const ycorsika = SVector{3}([0.89192975455881607, 0.18563051261662877, -0.41231374670066206])
const xcorsika = SVector{3}([0, -0.91184756344828699, -0.41052895273466672])
const zcorsika = whitepaper_normal_vec.proj
const xyzcorsika = inv([
    xcorsika.x xcorsika.y xcorsika.z;
    ycorsika.x ycorsika.y ycorsika.z;
    zcorsika.x zcorsika.y zcorsika.z;
])
const pavel_sim = jldopen("/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Oct16th2023_WhitePaper_300k.jld2")
const config = SimulationConfig(; pavel_sim["config"]...)
const injector = Tambo.Injector(config)
const geo = Tambo.Geometry(config)
const plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)

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
        "--first_n"
            help = "Only look at the first n events"
            arg_type = Int
            default = 0
    end
    return parse_args(s)
end

function add_hits!(d::Dict, events, modules)
    for event in events
        a = inside.(modules, Ref(event.pos))
        s = sum(a)
        @assert s <= 1 "Seems like you're in more than one module.. That doesn't seem right"
        if s > 0
            mod = modules[findfirst(a)]
            if !(mod.idx in keys(d))
                d[mod.idx] = 0
            end
            weight = event.weight
            if weight != 1
                weight = rand(Poisson(weight))
            end
            d[mod.idx] += weight
        end
    end
end

function has_triggered(hits_dict::Dict)
    d_filter = filter(x->x[2] >= 3, hits_dict) # Filter out modules below threshhold
    if length(d_filter) < 3
        return false
    end
    return sum(values(d_filter)) >= 30
end

function trigger_function(
    event_number,
    modules,
    basedir;
    xmax=nothing,
    xmin=nothing,
    ymin=nothing,
    ymax=nothing
)
    if isnothing(xmin)
        xmin = minimum([m.x for m in modules])
    end
    if isnothing(xmax)
        xmax = maximum([m.x for m in modules])
    end
    if isnothing(ymin)
        ymin = minimum([m.y for m in modules])
    end
    if isnothing(ymax)
        ymax = maximum([m.y for m in modules])
    end

    d = Dict{Int, Int}()
    triggered = false
    files = find_extant_files(event_number, basedir)
    while ~triggered && length(files) > 0
        file = pop!(files)
        # This nastiness needs to be here cuz sometimes the file is empty(?)
        # and it will error in an annoying way. I think the `pqf = nothing` is superfluous
        # but I put it there before and I'm worried it will cause issues
        pqf = nothing
        try
            pqf = pq.ParquetFile(file)
        catch
            continue
        end
        for batch in pqf.iter_batches()
            events = loadcorsika(batch)
            events = filter(
                e -> xmin < e.pos.x && e.pos.x < xmax && ymin < e.pos.y && e.pos.y < ymax,
                events
            )
            add_hits!(d, events, modules)
            triggered = has_triggered(d)
            if triggered
                break
            end
            # These two lines avoid runaway RAM usage
            events = nothing
            GC.gc()
        end
    end
    return triggered
end

function main()
    args = parse_commandline()
    
    modules = Tambo.make_trianglearray(
        -2000units.m,
        3000units.m,
        -args["length"] * units.m/2,
        args["length"] * units.m/2,
        args["deltas"] * units.m,
        ϕ=whitepaper_normal_vec.ϕ
    )

    modules = filter(
        (m,) -> zmin < Tambo.plane_z(m.x, m.y, plane) < zmax, modules
    )

    #xmax = maximum([m.x for m in modules])
    #xmin = minimum([m.x for m in modules])
    #ymax = maximum([m.y for m in modules])
    #ymin = minimum([m.y for m in modules])

    trigger_events = Int[]
    event_numbers = get_event_numbers(args["basedir"])
    if args["first_n"] > 0
        event_numbers = event_numbers[1:args["first_n"]]
    end
    for (idx, event_number) in enumerate(event_numbers)
        @show event_number
        triggered = trigger_function(event_number, modules, args["basedir"])
        if triggered
            push!(trigger_events, event_number)
            @show trigger_events
        end
        if idx % 100 == 0
            @show length(trigger_events) / idx
        end
    end
    np.save(args["outfile"], trigger_events)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
