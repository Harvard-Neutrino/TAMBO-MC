module PowerLaws

export PowerLaw, sample

struct PowerLaw
    γ::Float64
    emin::Float64
    emax::Float64
    norm::Float64
    function PowerLaw(γ, emin; emax=Inf)
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

function sample(pl::PowerLaw)
    u = rand()
    if pl.γ==1
        val = pl.emax^u/(pl.emin^(u-1))
    else
        mg = 1-pl.γ
        val = (u*pl.emax^mg + (1-u)*pl.emin^mg)^(1/mg)
    end
    val
end

end