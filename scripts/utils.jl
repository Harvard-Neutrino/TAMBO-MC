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

