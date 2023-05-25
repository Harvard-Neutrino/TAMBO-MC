This is an example of how some of the major pieces of the TAMBO MC work, using Pluto.jl. In order to install the necessary packages for running this, start a julia REPL session and run:

```julia-repl
julia> using Pkg; Pkg.instantiate(".")

julia> using Pluto

julia> Pluto.run()
```

This will take you to a browser wherer you can load the notebooks by entering the path to them. I am trying to come up with a cleaner way of doing this, but this is what I got for now.

You can also heck out a static html version of these notebooks here:
[Geometries demo](./resources/Geometries_demo.jl.html)
