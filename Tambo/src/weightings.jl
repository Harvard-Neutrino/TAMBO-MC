function p_mc(
    event::InjectionEvent,
    pl::PowerLaw,
    xs::CrossSection,
    anglesampler::UniformAngularSampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry
)
    p = (
        probability(anglesampler) * 
        probability(injectionvolume) * 
        density(event.final_state.position, geo) / event.genX * 
        probability(xs, event.entry_state.energy, event.final_state.energy) * 
        probability(pl, event.initial_state.energy)
    )
    return p
end

function p_phys(
    event::InjectionEvent,
    xs::CrossSection,
    geo::Geometry
)
    p = (
        event.physX * units.Na / units.Miso *
        density(event.final_state.position, geo) / event.physX * 
        xs.diff_xs(event.entry_state.energy, event.final_state.energy)
    )
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