using .Samplers
#using QuadGK

function p_mc_energy(Eν::Float64, γ::Float64, Emin::Float64, Emax::Float64)
    
    pl = PowerLaw(γ,Emin,Emax)
    
    power_law(E) = pl(E)
    result, error = quadgk(power_law,Emin,Emax)
          
    return pl(Eν)/result

end

function p_mc_energy(event::InjectionEvent,config::Dict)
    return p_mc_energy(event.energy,config[:γ],config[:emin],config[:emax])
end

function p_mc_Ω(θmin::Float64,ϕmin::Float64,θmax::Float64,ϕmax::Float64)
    #θ is the polar angle, ϕ is the azimuthal angle 
    polar(x) = sin(x)
    
    #return quadgk(polar,θmin,θmax)[1] * (ϕmax-ϕmin)
    return 0 
end

function p_mc_Ω(config::Dict)
   return  p_mc_Ω(config[:θmin],config[:ϕmin],config[:θmax],config[:ϕmax])
end

function p_mc_localρ(x::Float64,y::Float64,z::Float64,geo::Geometry)
  
  if inside(x,y,z,geo) == true
      return units.ρrock0 
  else 
      return units.ρair0
  end

end

function p_mc_localρ(event::InjectionEvent,geo::Geometry)
    return p_mc_localρ(event.initial_state.position.data...,geo)
end

function p_mc_area(r::Float64)
    return pi*(r/units.meters)^2
end

function p_mc_area(config::Dict)
    return p_mc_area(config[:r_injection])
end

function p_mc_xsection()
  return 0 
end

