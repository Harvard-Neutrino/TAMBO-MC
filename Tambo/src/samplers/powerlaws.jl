struct PowerLaw
    γ::Float64
    emin::Float64
    emax::Float64
    norm::Float64
    function PowerLaw(γ, emin, emax, norm)
        """
        Returns a normalized PowerLaw
        """
        if γ >= 1
            @assert !isinf(emax)
        end
        return new(γ, emin, emax, norm)
    end
end

function PowerLaw(γ, emin, emax)
    if γ == 1
        norm = 1 / log(emax / emin)
    else
        γ > 1
        mg = 1 - γ
        norm = mg / (emax^mg - emin^mg)
    end
    return PowerLaw(γ, emin, emax, norm)
end

function Base.show(io::IO, pl::PowerLaw)
    print(
        io,
        "PowerLaw(γ=$(pl.γ), emin=$(pl.emin / units.GeV) GeV, emax=$(pl.emax / units.GeV) GeV)",
    )
end

function (pl::PowerLaw)(e)
    return pl.norm * e^(-pl.γ)
end

function Base.rand(pl::PowerLaw)
    u = Base.rand()
    if pl.γ == 1
        b = (pl.emax) .^ u
        a = (pl.emin) .^ (u .- 1)
        return b ./ a
    end
    mg = 1 - pl.γ
    val = (u .* pl.emax .^ mg + (1 .- u) .* pl.emin^mg) .^ (1 / mg)
    return val
end

function Base.rand(n::Int, pl::PowerLaw)
    u = Base.rand(n)
    if pl.γ == 1
        b = (pl.emax) .^ u
        a = (pl.emin) .^ (u .- 1)
        val = b ./ a
        return val
    end
    mg = 1 - pl.γ
    val = (u .* pl.emax .^ mg + (1 .- u) .* pl.emin^mg) .^ (1 / mg)
    return val
end
