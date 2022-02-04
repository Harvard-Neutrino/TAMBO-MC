module EventDisplays

push!(LOAD_PATH, @__DIR__)
import Position:TPoint
import Tracks
import Particles

#=
Object for tracking a track as it moves from start to end
useful for making gifs and such
=#
mutable struct Trayectory
    track::Track
    # Parameter by which to increment the track
    dλ::Float64
    # Vectors for holding the position history of track
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    current_λ
    function Trayectory(track::Track)
        new(track, 0.01, [], [], [], 0.0)
    end
    function Trayectory(track::Track, dλ::Float64)
        new(track, dλ, [], [], [], 0.0)
    end
end

function step!(t::Trayectory)
    # Get current position
    p = find_position(t.track, t.current_λ)
    # Add positions to the Trayectory history
    append!(t.x, p.x)
    append!(t.y, p.y)
    append!(t.z, p.z)
    # Increment the Trayectory forward
    t.current_λ += t.dλ
end

function split_for_plot(x, y, z)
    region_change = diff([is_in_mountain(x...) for x in zip(x, y, z)]).!=0
    split_i       = [i for i in 1:length(region_change) if region_change[i]==1]
    split_track   = []
    old_idx       = 1
    current_color = ifelse(is_in_mountain(x[1], y[1], z[1]), :red, :black)
    for new_idx in split_i
        push!(split_track, (x[old_idx:new_idx], y[old_idx:new_idx], z[old_idx:new_idx], current_color))
        old_idx = new_idx
        current_color = ifelse(current_color==:black, :red, :black)
    end
    push!(split_track, tuple(x[old_idx:end], y[old_idx:end], z[old_idx:end], current_color))
    split_track
end

function split_for_plot(t::Trayectory)
    split_for_plot(t.x, t.y, t.z)
end


end