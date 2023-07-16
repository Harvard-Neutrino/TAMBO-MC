function p_mc(
    event::InjectionEvent,
    pl::PowerLaw,
    xs::CrossSection,
    anglesampler::UniformAngularSampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry
)
    p = probability(anglesampler)
    p *= probability(injectionvolume)
    p *= density(event.final_state.position, geo) / event.genX
    p *= probability(xs, event.entry_state.energy, event.final_state.energy)
    p *= probability(pl, event.initial_state.energy)
    #p = (
    #    probability(anglesampler) * 
    #    probability(injectionvolume) * 
    #    density(event.final_state.position, geo) / event.genX * 
    #    probability(xs, event.entry_state.energy, event.final_state.energy) * 
    #    probability(pl, event.initial_state.energy)
    #)
    return p
end


function p_phys(
    event::InjectionEvent,
    xs::CrossSection,
    geo::Geometry
)
    Miso = (938.27208816units.MeV + 939.5654133units.MeV) / 2
    Na = 6.02214076e23
    p = (event.physX +event.genX) * Na / Miso
    p *= density(event.final_state.position, geo) / (event.physX + event.genX)
    p *= xs.differential_xs(event.entry_state.energy, event.final_state.energy)
    #p = (
    #    event.physX * Na / Miso *
    #    density(event.final_state.position, geo) / event.physX * 
    #    xs.differential_xs(event.entry_state.energy, event.final_state.energy)
    #)
    return p
end

function oneweight(
    event::InjectionEvent,
    xsgen::CrossSection,
    xsphys::CrossSection,
    pl::PowerLaw,
    anglesampler::UniformAngularSampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo
)
    w = p_phys(event, xsphys, geo) / p_mc(event, pl, xsgen, anglesampler, injectionvolume, geo)
    return w
end
