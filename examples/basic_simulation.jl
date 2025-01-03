import Pkg
Pkg.develop(path="../Tambo/")

using Tambo
using JLD2

sim = Simulation("../resources/default_config.toml")
sim.config["steering"]["nevent"] = 10000

inject_ν!(sim)
propagate_τ!(sim)

jldopen("./data/basic_simulation_output.jld2", "w") do jldf
    jldf["simulation"] = sim
end
