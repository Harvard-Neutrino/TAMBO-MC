using Distributions: Uniform

abstract type AbstractInjectionShape end

include("crosssections.jl")
include("angularsamplers.jl")
include("injectionvolumes.jl")
include("powerlaws.jl")
include("injectionplane.jl")
