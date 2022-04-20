module Particles

using PowerLaws: PowerLaw, sample

struct Particle
    pdg_id::Int
    energy::Float64
    function Particle(pdg_id::Int, inv_cdf::Function)
        u = rand()
        e = inv_cdf(u)
        new(pdg_id, e)
    end
    function Particle(pdg_id::Int, pl::PowerLaw)
        # TODO this is sloppy
        e = sample(pl)[1]
        new(pdg_id, e)
    end
end

function interact(p::Particle)
    nothing
end

end