using Glob

using FileWatching: watch_file
using Base.Filesystem

function write_to_shared_file(sharedfilename, key, contenttowrite)
    lockacquired = false
    lockfilename = sharedfilename * ".lock"
    #lockfilename = lockfile_location * basename(sharedfilename) * ".lock"
    local lockfilehandle
    while !lockacquired
        while isfile(lockfilename)
            # watch_file will notify if the file status changes, waiting until then
            # here we want to wait for the file to get deleted
            watch_file(lockfilename)
        end
        try
            # try to acquire the lock by creating lock file with JL_O_EXCL (exclusive)
            lockfilehandle = Filesystem.open(lockfilename, JL_O_CREAT | JL_O_EXCL, 0o600)
            lockacquired = true
        catch err
            # in case the file was created between our `isfile` check above and the
            # `Filesystem.open` call, we'll get an IOError with error code `UV_EEXIST`.
            # In that case, we loop and try again.
            if err isa IOError && err.code == Base.UV_EEXIST
                continue
            else
                rethrow()
            end
        end
    end
    # now that the lock is acquired, we can safely write to
    # the actual shared file we want to write to
    jldopen(sharedfilename, "r+") do jldf
    #open(sharedfilename, append = true) do sharedfile
        # write to the sharedfile here just as usual
        jldf[key] = contenttowrite
    end

    # free up the lock so that the other process can acquire it if it needs
    close(lockfilehandle)
    Filesystem.unlink(lockfilename)
end

function get_event_numbers(basedir::String, simset_id::String)
    sims = glob("showers/$(simset_id)/shower_*_?", basedir)
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

function find_extant_files(simset_id::String, event_number::Int, basedir::String) Vector{String}
    files = glob("showers/$(simset_id)/shower_$(event_number)_?/*/particles.parquet", basedir)
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
    
    dir = (Matrix(df[:,["nx","ny","nz"]]) * xyzcorsika') * RotZ(-π/2)
    df.nx = dir[:,1] # particle direction in x
    df.ny = dir[:,2] # particle direction in y
    df.nz = dir[:,3] # particle direction in z

    ts = df.time .* units.second # time
    ws = df.weight
    return df
end 

"""
    particle_plane_intersection(p0::SVector{3}, dir::SVector{3}, pdet::SVector{3}, ndet::SVector{3})

Calculate the intersection point of a particle trajectory with a detector plane.

This function computes where a particle's trajectory (defined by its position and direction upon impacting the CORSIKA plane)
intersects with the plane of a detector (defined by a point on the plane and its normal vector).

# Arguments
- `p0::SVector{3}`: Position on the particle trajectory (3D vector)
- `dir::SVector{3}`: Direction vector of the particle trajectory (should be normalized)
- `pdet::SVector{3}`: A point on the detector plane (3D vector)
- `ndet::SVector{3}`: Normal vector to the detector plane (should be normalized)

# Returns
- `SVector{3}`: The 3D intersection point if the trajectory intersects the plane
- `nothing`: If the trajectory is parallel to the plane (no intersection)

# Mathematical Details
The intersection is calculated using the parametric line equation and plane equation:
- Line: `P(t) = p0 + t * dir`
- Plane: `dot(P - pdet, ndet) = 0`

Solving for parameter `t` gives: `t = dot(pdet - p0, ndet) / dot(dir, ndet)`
"""
function particle_plane_intersection(p0::SVector{3}, dir::SVector{3}, pdet::SVector{3}, nhat_det::SVector{3})
           denom = dot(dir, nhat_det)
           if abs(denom) < 1e-10
               return nothing  # Parallel, no intersection
           end
           t = dot(pdet - p0, nhat_det) / denom
           # TODO: check if t is "big." If so, particle may not have yet been created by shower
           return p0 + t * dir
       end

function detector_local_coords(p, center, nhat_det, uaxis)
    u = normalize(uaxis)
    v = normalize(cross(nhat_det, u))
    rel = p - center
    return (dot(rel, u), dot(rel, v))
end

function particle_hits_detector(
    origin::SVector{3}, direction::SVector{3},
    center::SVector{3}, n::SVector{3},
    width::Float64, height::Float64, uaxis::SVector{3}
)
    # 1. Intersection with detector plane
    hitpt = particle_plane_intersection(origin, direction, center, n)
    if hitpt === nothing
        return false
    end
    #println("Hit point: ", hitpt)
 
    # 2. Detector local coordinates
    (u, v) = detector_local_coords(hitpt, center, n, uaxis)
    return abs(u) <= width/2 && abs(v) <= height/2
end