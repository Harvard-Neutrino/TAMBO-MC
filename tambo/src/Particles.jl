module Particles

struct Particle
    pdg_id::Int
    energy::Float64
end

function interact(p::Particle)
    nothing
end

end