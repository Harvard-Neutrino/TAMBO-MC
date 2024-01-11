using Pkg
Pkg.activate("../Tambo")
using Tambo
using PyCall
using JLD2
using StaticArrays
using ArgParse
using Glob
using Distributions

ak = pyimport("awkward")
np = pyimport("numpy")

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

function loadcorsika(files)
    events = Tambo.CorsikaEvent[]
    for file in files
        x = nothing
        try
            x = ak.from_parquet(file)
        catch
            continue
        end
        xs = x["x"].to_numpy() .* units.m
        ys = x["y"].to_numpy() .* units.m
        zs = x["z"].to_numpy() .* units.m
        new_poss = []
        for (x, y, z) in zip(xs, ys, zs)
            push!(new_poss, SVector{3}(xyzcorsika * [x,y,z]))
        end
        ts = x["time"].to_numpy() .* units.second
        ws = x["weight"].to_numpy() .* 1.0
        ids = x["pdg"].to_numpy()
        es = x["kinetic_energy"].to_numpy() .* units.GeV
        for (id, e, pos, t, w) in zip(ids, es, new_poss, ts, ws)
            push!(events, Tambo.CorsikaEvent(id, e, pos, t, w))
        end
    end
    return events
end

function make_hits_dict(events, modules)
    d = Dict{Int, Int}()
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
    return d
end

function did_trigger_fast(events, modules)
    d = Dict{Int, Int}()
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
            if trigger_helper(d)
                return true
            end
        end
    end
    return false
end

function trigger_helper(d)
    d_filter = filter(x->x[2] >= 3, d) # Filter out modules below threshhold
    if length(d_filter) < 3
        return false
    end
    return sum(values(d_filter)) > 30
end

function did_trigger(events, modules)
    d = make_hits_dict(events, modules)
    return trigger_helper(d)
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

    println(length(modules))
    trigger_events = Int[]
    event_numbers = get_event_numbers(args["basedir"])
    for (idx, event_number) in enumerate(event_numbers)
        #@show event_number
        files = find_extant_files(event_number, args["basedir"])
        #@show files
        events = loadcorsika(files)
        #@show length(events)
        if did_trigger_fast(events, modules)
        #if did_trigger(events, modules)
            push!(trigger_events, event_number)
            println(length(trigger_events) / idx)
        end
        if idx % 100 == 0
            @show length(trigger_events) / idx
        end
    end
    np.save(args["outfile"], trigger_events)
end

main()
