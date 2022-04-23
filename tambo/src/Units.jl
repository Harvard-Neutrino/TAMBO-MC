module Units

using NaturallyUnitful

export GeV, MeV, m, cm, km, ρ_air0, ρ_rock0, mwe

GeV = u"GeV"
MeV = u"MeV"
PeV = u"PeV"
m = u"m"
cm = u"cm"
km = u"km"
mwe = 100u"g"/u"cm"^2
ρ_air0 = 1.225e-3u"g"/u"cm"^3
ρ_rock0 = 2.6u"g"/u"cm"^3

end # module