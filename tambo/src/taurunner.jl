using PyCall
tr = pyimport("taurunner")
pp = pyimport("proposal")
np = pyimport("numpy")

nevents   = 1000
eini      = 1e19  #initial energy (eV)
theta     = 89.0  #degrees (nadir) 
pid       = 16    #flavor
xs_model  = "CSMS"

xs        = tr.cross_sections.CrossSections(xs_model)
#earth     = tr.body.Body(1.0, 1e6)
earth     = tr.body.earth.construct_earth(layers=[(4., 1.0)])
thetas    = tr.utils.make_initial_thetas(nevents, theta)
energies  = tr.utils.make_initial_e(nevents, eini)
prim_prop = tr.utils.make_propagator(pid, earth)
rand      = np.random.RandomState(seed=8)

output = tr.main.run_MC(energies, 
                thetas, 
                earth, 
                xs,
                prim_prop, 
                rand,)

struct TROutput
    Eini::Array{Float64}
    Eout::Array{Float64}
    Theta::Array{Float64}
    nCC::Array{Int}
    nNC::Array{Int}
    PDG_Encoding::Array{Int}
    event_ID::Array{Int}
    final_position::Array{Float64}
end

Eini = output.__getitem__("Eini")
Eout = output.__getitem__("Eout")
Theta = output.__getitem__("Theta")
nCC = output.__getitem__("nCC")
nNC = output.__getitem__("nNC")
PDG_Encoding = output.__getitem__("PDG_Encoding")
event_ID = output.__getitem__("event_ID")
final_position = output.__getitem__("final_position")
out = TROutput(Eini, Eout, Theta, nCC, nNC, PDG_Encoding, event_ID, final_position)

mask = out.PDG_Encoding.==15

#propagate taus in air
atmosphere = tr.body.Body(0.001225, 100000.)   #body(density [g/cm^3], radius [km])
air_prop   = tr.utils.make_propagator(pid, atmosphere)
proposal_lep           = pp.particle.DynamicData(pp.particle.TauMinusDef().particle_type)
proposal_lep.position  = pp.Vector3D(0, 0, 0)
proposal_lep.direction = pp.Vector3D(0, 0, 1)
track  = tr.track.Chord(theta=0., depth=0.)

tau_out = []
for tau in exiting_taus
    this_tau = tr.particle.Particle(tau['PDG_Encoding'], tau['Eout'], 0., 0., 0, xs, air_prop, proposal_lep, False, False)
    this_tau.PropagateChargedLepton(atmosphere, track)
    #energy at decay in eV, distance to decay in km
    tau_out.append((this_tau.energy, this_tau.chargedposition))
end
