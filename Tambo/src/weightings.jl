function p_mc(
    event::InjectionEvent,
    pl::PowerLaw,
    xs::CrossSection,
    anglesampler::UniformAngularSampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry
)
    p = probability(anglesampler, event)
    p *= probability(injectionvolume, event)
    p *= density(event.final_state.position, geo) / event.genX
    p *= probability(xs, event)
    p *= probability(pl, event)
    return p
end

function p_phys(
    event::InjectionEvent,
    xs::CrossSection,
    geo::Geometry
)
    Miso = (938.27208816units.MeV + 939.5654133units.MeV) / 2
    p = event.genX / Miso
    p *= density(event.final_state.position, geo) / event.genX
    p *= xs.differential_xs(event.entry_state.energy, event.final_state.energy)
    return p
end

function oneweight(
    event::InjectionEvent,
    xsgen::CrossSection,
    xsphys::CrossSection,
    pl::PowerLaw,
    anglesampler::UniformAngularSampler,
    injectionvolume::SymmetricInjectionCylinder,
    geo::Geometry
)
    n = p_phys(event, xsphys, geo)
    d = p_mc(event, pl, xsgen, anglesampler, injectionvolume, geo)
    if n==0
        return 0
    end
    w = n / d
    return w
end

function oneweight(
    event::InjectionEvent,
    injector::Injector,
    xsphys::CrossSection,
)
    return oneweight(event, injector.xs, xsphys, injector.powerlaw, injector.anglesampler, injector.injectionvolume, injector.geo)
end
