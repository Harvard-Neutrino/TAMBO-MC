module CrossSections

export outgoing_energy

using PyCall
push!(LOAD_PATH, @__DIR__)
using Units: eV

py"""
import pickle
 
def load_pickle(fpath):
    with open(fpath, "rb") as f:
        data = pickle.load(f)
    return data
"""

#tr = pyimport("taurunner")
#xs_path = tr.__path__ * "/resources/cross_section_tables/"
xs_path = "/Users/jlazar/research/TauRunner/taurunner/resources/cross_section_tables/"
println(xs_path)

load_pickle = py"load_pickle"
nu_diff_xs = load_pickle(xs_path * "/CSMS_nu_p_dsde_CC.pkl")
nubar_diff_xs = load_pickle(xs_path * "CSMS_nubar_p_dsde_CC.pkl")
diff_xs = [nu_diff_xs, nubar_diff_xs]
EMIN = 1e9
energy_fractions = LinRange(0, 1, 1000)[2:end-1]

function call_diff_xs(ein, zz, spl)
    ein_ev = (ein |> eV).val
    println(ein_ev)
    println((log(ein_ev), zz))
    #res = exp(spl(log(ein_ev), zz)[1])/ein_ev
end

function outgoing_energy(eν, νtype, n)
    w_int = ProbabilityWeights(call_diff_xs(eν, energy_fractions, diff_xs[νtype]))
    z_choice = StatsBase.sample(energy_fractions, w_int, n)
    println(z_choice)
    eτ = z_choice*(eν-EMIN)+EMIN
end

end # module