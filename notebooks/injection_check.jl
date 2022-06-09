### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ ec8a942c-cde8-11ec-2793-edb806aa17a0
begin
	using Pkg
	Pkg.activate("/Users/jlazar/research/TAMBO-MC/tambo/")
	using Tambo
	using HDF5
end

# ╔═╡ 977eb838-809b-4325-8d73-109f688b109f
begin
	using PlutoUI
	using Plots
	using Printf
	units = Tambo.units
end

# ╔═╡ e6f9573b-4aa0-4d71-b701-89f87b625ed8
begin
	fname = "/Users/jlazar/research/TAMBO-MC/10000_new_injection_check.h5"
	fid = h5open(fname, "r");
end

# ╔═╡ 4b0bcae6-dd2e-4d6a-8747-4fc56e668ec9
begin
	ts = TAMBOSim();
	ts.n = 100
	# ts()
end

# ╔═╡ a8c31968-0d6a-4259-8587-d33bd2970420
@bind k Slider(1:50; show_value=true)

# ╔═╡ 202b11c2-90ca-4306-b49d-de2b90fbc25b
begin
	xx = fid["interaction_vertex_x"][]
	yy = fid["interaction_vertex_y"][]
	zz = fid["interaction_vertex_z"][]
	xrange = LinRange(minimum(xx), maximum(xx), 51)
	yrange = LinRange(-20000, 20000, 200)
	xl, xh = xrange[k], xrange[k+1]
	mask = xx .<= xh .&& xx .>= xl
	plt = scatter(yy[mask], zz[mask], alpha=0.1, label=false)
	plot!(plt, xlim=(minimum(yy), maximum(yy)), ylim=(minimum(zz), 5000))
    plot!(plt, yrange, ts.geo.(Ref((xl+xh)*units[:m]/2), yrange.*units[:m])./units[:m], label="mid")
    plot!(plt, yrange, ts.geo.(Ref(xh*units[:m]), yrange.*units[:m])./units[:m], label="high")
    plot!(plt, yrange, ts.geo.(Ref(xl*units[:m]), yrange.*units[:m])./units[:m], label="low")
	yrange = LinRange(-20000, 20000, 200)
	plot!(plt, xlabel="Y [m]", ylabel="Z [m]")
    string_xl = @sprintf("%5.0f",xl)
    string_xh = @sprintf("%5.0f",xh)
    string_xlxh = string_xl * " m<=x<" * string_xh*" m"
    annotate!(plt, -1.7e4, -2.75e4, text(string_xlxh, :left, 9, "sans_serif"))
end

# ╔═╡ c888893d-107c-4389-a2de-eb704f304dae


# ╔═╡ 5b5b8c76-fcbc-40ec-af76-6a610b712aea
@bind j Slider(0:99, show_value=true)

# ╔═╡ ae29d253-babb-48fa-8cb6-0c1f88850c84
@bind i Slider(0:100, show_value=true)

# ╔═╡ 356d29e1-dbb6-499b-b00f-12e6b946d72d
begin
	pnxx = fid["p_near_x"][]
	pnyy = fid["p_near_y"][]
	pnzz = fid["p_near_z"][]
	θ = π * j / 99
	ϕ = 2π * i / 100
	xn = sin(θ)cos(ϕ)
	yn = sin(θ)sin(ϕ)
	zn = cos(θ)
	# See if Julia has a way to plot planes
	z(x, y) = (x * xn + y * yn) / zn
	pnxrange = pnyrange = LinRange(-900, 900, 100)
	# l = @layout [a{0.5w}]
	plt1 = surface(pnxrange, pnyrange, z, color="black", alpha=0.0, colorbar=false)
	plot!(plt1, xlim=(-900, 900), ylim=(-900, 900), zlim=(-900, 900))
	scatter!(plt1, pnxx, pnyy, pnzz, color="blue", alpha=0.1)
	plot!(plt1, xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")
	plot!(camera=(30,30))
end

# ╔═╡ aac01ef6-67fb-4978-a49c-b831ec879e5e


# ╔═╡ 2cd87d15-5d4a-40de-ade8-e9d44cdc9c15
ts.geo.box.c2./units[:m]

# ╔═╡ Cell order:
# ╠═ec8a942c-cde8-11ec-2793-edb806aa17a0
# ╠═e6f9573b-4aa0-4d71-b701-89f87b625ed8
# ╠═4b0bcae6-dd2e-4d6a-8747-4fc56e668ec9
# ╠═977eb838-809b-4325-8d73-109f688b109f
# ╟─a8c31968-0d6a-4259-8587-d33bd2970420
# ╟─202b11c2-90ca-4306-b49d-de2b90fbc25b
# ╠═c888893d-107c-4389-a2de-eb704f304dae
# ╟─5b5b8c76-fcbc-40ec-af76-6a610b712aea
# ╟─ae29d253-babb-48fa-8cb6-0c1f88850c84
# ╠═356d29e1-dbb6-499b-b00f-12e6b946d72d
# ╠═aac01ef6-67fb-4978-a49c-b831ec879e5e
# ╠═2cd87d15-5d4a-40de-ade8-e9d44cdc9c15
