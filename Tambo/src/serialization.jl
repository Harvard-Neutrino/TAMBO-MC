function flatten_obj(obj; pre="", lists=nothing)
    if typeof(lists)==Nothing
        lists = Tuple{String, Int}[]
    end
    for field in fieldnames(typeof(obj))
        val = getfield(obj, field)
        t = typeof(val)
        if t <: Real
            s = length(pre)==0 ? string(field) : "$(pre).$(string(field))"
            push!(lists, (s, 0))
            continue
        elseif val isa AbstractVector
            s = length(pre)==0 ? string(field) : "$(pre).$(string(field))"
            push!(lists, (s, length(val)))
            continue
        else
            w = length(pre)==0 ? string(field) : "$(pre).$(string(field))"
            flatten_dict(val; pre=w, lists=lists)
        end
    end
    return lists
end

function write_to_h5(objs, baseg)
    types = flatten_obj(first(objs))

    for (s, l) in types
        out = l==0 ? zeros(length(objs)) : zeros((length(objs), l))
        for (idx, event) in enumerate(objs)
            v = recursively_access(event, s)
            if l==0
                out[idx] = v
            else
                out[idx, :] = v
            end
        end

        splits = split(s, ".")
        g = baseg
        for x in splits[1:end-1]
            if ~(x ∈ keys(g))
                create_group(g, x)
            end
            g = g[x]
        end
        create_dataset(g, splits[end], out)
    end
end

function recursively_access(obj, s)
    for x in split(s, ".")
        obj = getfield(obj, Symbol(x))
    end
    return obj
end
