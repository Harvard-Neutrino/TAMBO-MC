using Distributions
using Random

function accept_all(e) :: Float64
    return 1.0
end

function get_nhit(hit_map, efficiencies)
#function get_nhit(hit_map, efficiencies::Dict{Int, Function})
    nhit = 0
    for hits in values(hit_map)
        nhit += sum([corsika_int_weight(hit, efficiencies) for hit in hits])
    end
    return nhit
end

function did_trigger(
    hit_map,
    module_trigger_thresh::Int=3,
    event_trigger_thresh::Int=30
)
    pids = Int[]
    for hits in values(hit_map)
        for hit in hits
            if hit.pdg in pids
                continue
            end
            push!(pids, hit.pdg)
        end
    end
    efficiencies = Dict(pid => accept_all for pid in pids)

    return did_trigger(
        hit_map,
        efficiencies,
        module_trigger_thresh,
        event_trigger_thresh
    )
end

function did_trigger(
    hit_map,
    efficiencies,
    #efficiencies::Dict{Int, Function},
    module_trigger_thresh::Int=3,
    event_trigger_thresh::Int=30
)
    hit_map = filter(
        kv->sum(corsika_int_weight(kv[2], efficiencies)) >= module_trigger_thresh,
        hit_map
    )
    if length(hit_map) < 3
        return false
    end
    nhit = get_nhit(hit_map, efficiencies)
    return nhit >= event_trigger_thresh
end

function corsika_int_weight(
    event::Tambo.CorsikaEvent,
    efficiencies
)::Int
    seed = mod(hash(event), 2^32)
    return corsika_int_weight(event, efficiencies, seed)
end

function corsika_int_weight(
    event::Tambo.CorsikaEvent,
    efficiencies,
    seed::Int
) :: Int
    Random.seed!(seed)
    efficiency = efficiencies[event.pdg]
    ϵ = efficiency(event.kinetic_energy)
    if rand() > ϵ
        return 0
    end
    weight = event.weight
    if weight > 1
        weight = rand(Poisson(weight))
    end
    return weight
end

function corsika_int_weight(
    events::Vector{Tambo.CorsikaEvent},
    efficiencies
) :: Vector{Int}
    return [corsika_int_weight(event, efficiencies) for event in events]
end
