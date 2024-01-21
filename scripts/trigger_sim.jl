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
const pavel_sim = jldopen("/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/WhitePaper_300k.jld2")
const config = SimulationConfig(; pavel_sim["config"]...)
const injector = Tambo.Injector(config)
const geo = Tambo.Geometry(config)
const plane = Tambo.Plane(whitepaper_normal_vec, whitepaper_coord, geo)

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
    end
    return parse_args(s)
end

function loadcorsika(batch)
    ids = np.array(batch[2])
    es = np.array(batch[3]) .* units.GeV
    xs = np.array(batch[4]) .* units.m # x coord
    ys = np.array(batch[5]) .* units.m # y coord
    zs = np.array(batch[6]) .* units.m # z coord
    poss = []
    for (x, y, z) in zip(xs, ys, zs)
        push!(poss, SVector{3}(xyzcorsika * [x,y,z]))
    end
    ts = np.array(batch[10]) .* units.second # time
    ws = np.array(batch[11])
    events = Tambo.CorsikaEvent.(ids, es, poss, ts, ws)
    return events
end

function find_hits!(d::Dict, events, modules)
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

function did_trigger_fast(events, modules)
    d = Dict{Int, Int}()
    for event in events
        a = inside.(modules, Ref(event.pos))
        s = sum(a)
        @assert s <= 1 "Seems like you're in more than one module. That doesn't seem right"
        if s == 0
            continue
        end

        mod = modules[findfirst(a)]
        if !(mod.idx in keys(d))
            d[mod.idx] = 0
        end

        weight = event.weight
        if weight != 1
            weight = rand(Poisson(weight))
        end
        d[mod.idx] += weight
        if has_triggered(d)
            return true
        end
    end
    return false
end

function has_triggered(d)
    d_filter = filter(x->x[2] >= 3, d) # Filter out modules below threshhold
    if length(d_filter) < 3
        return false
    end
    return sum(values(d_filter)) > 30
end

function did_trigger(events, modules)
    d = make_hits_dict(events, modules)
    return has_triggered(d)
end

function get_event_numbers(basedir)
    sims = glob("sim_test_*_?", basedir)
    event_numbers = []
    for sim in sims
        event_number = parse(Int, split(sim, "_")[end-1])
        if event_number in event_numbers
            continue
        end
        push!(event_numbers, event_number)
    end
    return sort(event_numbers)
end

function find_extant_files(run_number::Int, basedir::String) Vector{String}
    particle_number = 1
    files = String[]
    while true # iterate until the file doesn't exist
        path = "$(basedir)/sim_test_$(run_number)_$(particle_number)/particles/particles.parquet"
        if !ispath(path)
          break
        end
        push!(files, path)
        particle_number += 1
    end
    return files
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
    xmax = maximum([m.x for m in modules])
    xmin = minimum([m.x for m in modules])
    ymax = maximum([m.y for m in modules])
    ymin = minimum([m.y for m in modules])

    trigger_events = Int[]
    event_numbers = get_event_numbers(args["basedir"])
    for (idx, event_number) in enumerate(event_numbers)
        @show event_number
        d = Dict{Int, Int}()
        triggered = false
        files = find_extant_files(event_number, args["basedir"])
        while ~triggered && length(files) > 0
            file = pop!(files)
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
                find_hits!(d, events, modules)
                events = nothing
                triggered = has_triggered(d)
                if triggered
                    break
                end
                GC.gc()
            end
        end
        if triggered
            push!(trigger_events, event_number)
            @show trigger_events
            #@show length(trigger_events) / idx
        end
        if idx % 100 == 0
            @show length(trigger_events) / idx
        end
    end
    np.save(args["outfile"], trigger_events)
end

main()
