module Samplers

using Distributions: Uniform

export UniformAngularSampler
export SymmetricInjectionCylinder, AsymmetricInjectionCylinder
export PowerLaw
export OutgoingCCEnergy

include("crosssections.jl")
include("angularsamplers.jl")
include("injectionvolumes.jl")
include("powerlaws.jl")

end # module
