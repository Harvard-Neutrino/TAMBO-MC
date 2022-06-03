This will introduce the major components of the TAMBO-MC. All REPL command are written assuming you launched it in this directory

# Geometry

The geometry contains two major `struct`s, `TPoint` and `Box`. The former is just a cartesian point that has some convenient functionalities added on to it. Idon't know if it is extremely necessary at this point, but we have it and it works so what can you do. The latter defines the generation region for the simulation. They can be used in the following manner.

```julia
using Pkg; Pkg.activate(".")
# Tell the REPL where our code lives
push!(LOAD_PATH, String(@__DIR__)*"/src")
# Load our module
using Geometry
#=
Create a box with corners at (30, 30, 140) and (60, 60, 100) and with
edges parallel to the principal axes. All boxes have edges parallel to the
principle axes !!
=#
b = Box((30,30,140), (60, 60, 100))
# We can sample points within trhe box
xyz = sample_xyz(b)
# And create a TPoint inside the box
p = TPoint(xyz)
# Check the it is inside the box
is_inside(p, b)
```

I have set aside space in the module for a `GenerationRegion` `struct` but there are a couple steps left before I can make it work

# Tracks

This module contains the `Track` and `Direction` `struct`s as well as some functions to act upon them. The former `struct` defines the trjectory of a particle, which can be constructed in a number of ways. For instance, you can give it a begining `TPoint` and an ending `TPoint`

```julia
using Pkg; Pkg.activate(".")
# Tell the REPL where our code lives
push!(LOAD_PATH, String(@__DIR__)*"/src")
# Load our module
using Tracks
ipoint = TPoint(15, 25, 50)
fpoint = TPoint(45, 55, 50)
# Make a track from these. Note this has no z-component
t = Track(ipoint, fpoint)
# Check that it is moving horizontally
t.direction.θ==π/2
```

We can also instantiate it using an initial point, a direction, and requiring the it ends on a supplied box

```julia
using Pkg; Pkg.activate(".")
# Tell the REPL where our code lives
push!(LOAD_PATH, String(@__DIR__)*"/src")
# Load our modules
using Tracks, Geometry
b = Box((30,30,140), (60, 60, 100))
xyz = sample_xyz(b)
ipoint = TPoint(xyz)
# Randomly sample a direction
d = Direction(rand()*π, rand()*2*π)
t = Track(ipoint, d, b)
# Check that it does in fact end on the box
t(1)
#=
We usually want tracks to enter from the box.
We can achieve this with the `reverse` function
=#
t = reverse(t)
# check that it starts on an edge of the box
t(0)
# Check that it ends at the initial point
t(1)==ipoint
```

This module also has a function for calculating the total column depth given a function. You can check this out in the `/notebooks/event_displays.ipynb`
