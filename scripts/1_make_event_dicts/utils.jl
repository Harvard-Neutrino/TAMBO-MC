using Glob

function get_event_numbers(basedir::String, simset::String, subsimset::String)
    sims = glob("showers/$(simset)_$(subsimset)/shower_*_?", basedir)
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

function find_extant_files(simset::String, subsimset::String, event_number::Int, basedir::String) Vector{String}
    files = glob("showers/$(simset)_$(subsimset)/shower_$(event_number)_?/*/particles.parquet", basedir)
    if isempty(files) 
        println("No files matched the pattern. The array is empty.")
    end 
    return files
end

function loadcorsika(df,xyzcorsika)
    ids = df.pdg
    es = df.kinetic_energy .* units.GeV
    pos = (Matrix(df[:,["x","y","z"]]) * xyzcorsika') * RotZ(-π/2)
    #pos = Matrix(df[:,["x","y","z"]]) * xyzcorsika'
    df.x = pos[:,1] .* units.m # x coord
    df.y = pos[:,2] .* units.m # y coord
    df.z = pos[:,3] .* units.m # z coord
    ts = df.time .* units.second # time
    ws = df.weight
    #events = Tambo.CorsikaEvent.(ids, es, pos, ts, ws)
    return df
end 