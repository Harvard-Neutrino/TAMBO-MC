using Pkg
Pkg.activate("../../Tambo")
using Tambo
using JLD2

include("trigger_defs.jl")

interpolated_effs = load("../../resources/detector_efficiencies/initial_IceTop_panel_interpolations.jld2")
global interpolated_eff_gamma = interpolated_effs["gamma_interp"]
global interpolated_eff_muon = interpolated_effs["muon_interp"]
global interpolated_eff_electron = interpolated_effs["electron_interp"]

trigger_type = "whitepaper"
module_trigger_thresh = 3
event_trigger_thresh = 30

#trigger_type = "threeamigos"
#module_trigger_thresh = 3
#event_trigger_thresh = 30

#trigger_type = "icetop_tanks"
#module_trigger_thresh = 300
#event_trigger_thresh = 3000

#trigger_type = "icetop_panels"
#module_trigger_thresh = 3 * 65
#event_trigger_thresh = 30 * 65

TAMBO_PATH = "../../Tambo"
# This should be the numpy file with the indicess of the triggered events
#EVENT_DICTS_PATH = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will/event_dicts/Jan7th2024_WhitePaper_300k_no_thin/"
EVENT_DICTS_PATH = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/will_misc/triggered_events/Jan7th2024_WhitePaper_300k_no_thin/"
# This should be the simulation parameters jld2 file
SIMULATION_FILE = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/corsika8/corsika-work/Oct16th2023_WhitePaper_300k.jld2"

sim = jldopen(SIMULATION_FILE)
config = SimulationConfig(; geo_spline_path="/n/home02/thomwg11/tambo/TAMBO-MC/resources/tambo_spline.jld2", filter(x->x[1]!=:geo_spline_path, sim["config"])...)
geo = Tambo.Geometry(config)
injector = Tambo.Injector(config)

event_dicts = load(EVENT_DICTS_PATH*"test_1_100.jld2")
for i in 2:100
    merge!(event_dicts, load(EVENT_DICTS_PATH*"test_$(i)_100.jld2"))
end

triggered_event_ids = []
for (key, value) in event_dicts
    if did_trigger(value, module_trigger_thresh, event_trigger_thresh, trigger_type)
        push!(triggered_event_ids, parse(Int, split(key, "/")[end]))
    end
end
triggered_events = sim["injected_events"][triggered_event_ids]

save("triggered_events_$(trigger_type)_mod_$(module_trigger_thresh)_event_$(event_trigger_thresh).jld2", "triggered_events", triggered_events)
