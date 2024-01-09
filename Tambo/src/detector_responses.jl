struct DetectionModule
  x::Float64
  y::Float64
  Δx::Float64
  Δy::Float64
  idx::Int
end

# This is not really true... Fix this to include events and not just indices
struct Hit
  mod::DetectionModule
  event::CorsikaEvent
end

function make_trianglearray(x0::Number, x1::Number, y0::Number, y1::Number, ds::Number; ϕ::Number=0.0)
  dy = ds * sqrt(3) / 2
  ys = y0:dy:y1
  dx = ds
  xs = x0:dx:x1
  m = [cos(ϕ) -sin(ϕ) 0 ; sin(ϕ) cos(ϕ) 0; 0 0 1]
  detection_modules = DetectionModule[]
  mod_idx = 1
  for x in xs for (idx, y) in enumerate(ys)
    offset = ds / 2
    if idx %2==0
       offset = 0
    end
    v = m * SVector{3}(x+offset, y, 0)
    mod = DetectionModule(v[1], v[2], 1, 1, mod_idx)
    push!(detection_modules, mod)
    mod_idx += 1
  end end
  return detection_modules
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

function compute_minimum_distance(xy::SVector{2}, detmods::Vector{DetectionModule})
  Δxs = xy.x .- [x.x for x in detmods]
  Δys = xy.y .- [x.y for x in detmods]
  x = sqrt.(Δxs .^2 .+ Δys .^2)
  idx = argmin(x)
  return idx, x[idx]
end

function compute_minimum_distance(event::CorsikaEvent, detmods::Vector{DetectionModule})
    xy = SVector{2}([event.x, event.y])
    return compute_minimum_distance(xy, detmods)
end

function plane_z(x::Number, y::Number, plane::Plane)
    z = (
         -plane.n̂.proj.x * (x - plane.x0.x)
         -plane.n̂.proj.y * (y - plane.x0.y)
    )
    z /= plane.n̂.proj.z
    return z
end

function find_near_hits(
  events::Vector{CorsikaEvent},
  detmods::Vector{DetectionModule};
  thresh=1units.m
)
  hits = Hit[]
  for (jdx, event) in enumerate(events)
    module_idx, dist = compute_minimum_distance(event, detmods)
    if dist > thresh
      continue
    end
    push!(
      hits,
      Hit(detmods[module_idx], event)
    )
  end
  return hits
end

function make_hit_dict(hits::Vector{Hit})
  d = Dict(hit.mod.idx => 0 for hit in hits)
  for hit in hits
    if hit.event.weight==1 # these events do not have thinning
      d[hit.mod.idx] += 1
    else # These do have thinning and thus we choose Poisson nunber
      d[hit.mod.idx] += rand(Poisson(hit.event.weight))
    end
  end
  return d
end
