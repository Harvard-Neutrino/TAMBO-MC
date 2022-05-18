module Units

using NaturallyUnitful

export GeV, MeV, m, cm, km, ρ_air0, ρ_rock0, mwe, eV, gr

GeV = 1e9
TeV = 1e12
PeV = 1e15
eV = 1
m = 5067730.93741
km = 1000.0 * m
cm = m / 100.0
gr = 5.62e+32
mwe = 100 * gr / cm^2
ρ_air0 = 1.225e-3 * gr / cm^3
ρ_rock0 = 2.6 * gr / cm^3

#GeV = u"GeV"
#MeV = u"MeV"
#PeV = u"PeV"
#eV = u"eV"
#m = u"m"
#cm = u"cm"
#km = u"km"
#mwe = 100u"g"/u"cm"^2
#ρ_air0 = 1.225e-3u"g"/u"cm"^3
#ρ_rock0 = 2.6u"g"/u"cm"^3

end # module