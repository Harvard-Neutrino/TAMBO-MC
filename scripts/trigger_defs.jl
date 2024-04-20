using Distributions

function get_nhit(hit_map)
    nhit = 0
    for hits in values(hit_map)
        nhit += sum([corsika_int_weight(hit) for hit in hits])
    end
    return nhit
end

function did_trigger(hit_map, module_trigger_thresh=3, event_trigger_thresh=30)
    if length(hit_map) < 3
        return false
    end
    hit_map = filter(kv->length(kv[2]) >= module_trigger_thresh, hit_map)
    if length(hit_map) < 3
        return false
    end
    nhit = get_nhit(hit_map)
    return nhit >= event_trigger_thresh
end

function corsika_int_weight(event::Tambo.CorsikaEvent)
    weight = event.weight
    if weight > 1
        weight = rand(Poisson(weight))
    end
    return weight
end