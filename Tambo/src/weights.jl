using .Samplers

function p_mc_energy(E_ν::Float64, γ::Float64, Emin::Float64, Emax::Float64)
    
    pl = PowerLaw(γ,Emin,Emax)
    return pl.norm      
end
