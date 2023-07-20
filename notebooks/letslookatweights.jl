### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 5c90a496-267b-11ee-39be-6b026a69e54b
begin 
	using Dierckx
	tambopath = realpath("/Users/pavelzhelnin/Documents/physics/TAMBO/tambo");
	using Pkg; 
	Pkg.activate(tambopath)
	
	using Tambo

end 

# ╔═╡ 9c37ad72-327a-470f-9c6e-04c861e745e8
begin 
	Pkg.add("PlutoUI")
	Pkg.add("Plots")
	using PlutoUI
	Pkg.add("CSV")
	Pkg.add("DataFrames")
	Pkg.add("LinearAlgebra")
	Pkg.add("Parquet2")
	using Plots,CSV,DataFrames,LinearAlgebra
	using StaticArrays
	using Statistics 
	using Parquet2:Dataset
	using Random
end

# ╔═╡ 03da697b-5efc-4562-9947-08844ba075bf
simulator = SimulationConfig("/Users/pavelzhelnin/Documents/physics/TAMBO/notebooks/fixed_file.jld2")

# ╔═╡ c5a2d246-a888-4555-a9a2-66749331b2d9
geo = Geometry("$(tambopath)/../resources/tambo_spline.jld2", simulator.tambo_coordinates);

# ╔═╡ b5058a30-243a-46a1-a8ff-faac41b8d7fa
injected_event = simulator["injected_events"]

# ╔═╡ 98de2ef7-dbf9-42cf-af59-a0f03cd2372a
Tambo.oneweight(injected_event[1])

# ╔═╡ Cell order:
# ╠═5c90a496-267b-11ee-39be-6b026a69e54b
# ╠═9c37ad72-327a-470f-9c6e-04c861e745e8
# ╠═03da697b-5efc-4562-9947-08844ba075bf
# ╠═c5a2d246-a888-4555-a9a2-66749331b2d9
# ╠═b5058a30-243a-46a1-a8ff-faac41b8d7fa
# ╠═98de2ef7-dbf9-42cf-af59-a0f03cd2372a
