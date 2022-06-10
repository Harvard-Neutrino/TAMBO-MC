struct PowerLaw
    γ::Float64
    emin::Float64
    emax::Float64
    norm
    function PowerLaw(γ, emin, emax)
        """
        Returns a normalized PowerLaw
        """
        if γ >= 1
            @assert !isinf(emax)
        end
        if γ==1
            norm = 1/log(emax/emin)
        else γ>1
            mg = 1-γ
            norm = mg/(emax^mg - emin^mg)
        end
        new(γ, emin, emax, norm)
    end
end

function (pl::PowerLaw)(e)
    pl.norm*e^(-pl.γ)
end

function Base.rand(pl::PowerLaw)
    u = Base.rand()
    if pl.γ==1
        b = (pl.emax) .^ u
        a = (pl.emin ) .^(u.-1)
        val = b ./ a 
    else
        mg = 1-pl.γ
        val = (u .* pl.emax .^ mg + (1 .- u) .* pl.emin^mg) .^ (1/mg)
    end
    val
end

function Base.rand(n::Int, pl::PowerLaw)
    u = Base.rand(n)
    if pl.γ==1
        b = (pl.emax) .^ u
        a = (pl.emin) .^(u.-1)
        val = b ./ a
    else
        mg = 1-pl.γ
        val = (u .* pl.emax .^ mg + (1 .- u) .* pl.emin^mg) .^ (1/mg)
    end
    val
end
