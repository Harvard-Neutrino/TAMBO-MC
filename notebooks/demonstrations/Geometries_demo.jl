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

# ╔═╡ d260f8b0-d780-11ec-1e84-695ebd7b2393
begin
	using Pkg
	Pkg.activate(".")
	tambo_path = "/Users/jlazar/research/TAMBO-MC/tambo/src/"
	push!(LOAD_PATH, tambo_path);
end

# ╔═╡ 78fb3deb-4202-4e91-96b3-b89fbd99eaf4
begin
	using StaticArrays
	using BenchmarkTools
end

# ╔═╡ 4b006b1e-f0bc-4afc-8544-4b1d44d29afe
begin
	using Geometries: Direction
	θ = π/3
	ϕ = π/4
	d = Direction(θ, ϕ)
end

# ╔═╡ 0312f8e0-50b5-41d0-9342-1d0835c4f846
begin
	using Geometries: Box
	c0 = [2,5,9]
	c1 = [-2, -5, -9]
	b = Box(c0, c1)
end

# ╔═╡ f28c2ca0-e241-4e02-bfc2-445688c81bc6
begin
	using Geometries: is_inside
	pt0 = SVector{3}([0, 0, 0])
	is_inside(pt0, b)
end

# ╔═╡ 38ac33c5-a9e6-40cd-9244-afd6547f6488
begin
	using Geometries: Geometry, load_spline
	spl = load_spline(tambo_path * "../../resources/tambo_spline.npy")
	geo1 = Geometry(spl)
end

# ╔═╡ c04d426a-0e52-4d5b-a0a2-0b16c4688e93
begin
	using Units
	geo1(0m, 0m) / m
end

# ╔═╡ 9ea8bd84-efa1-4f54-a077-bf8faa9698d4
using PlutoUI: Slider

# ╔═╡ 7ded89b6-496e-4bed-b3a2-53e8be42f6a3
md"""
Most of the TAMBO geometries rest on top of the `StaticArrays.SVector`. This is a type of vector where the length and type are fixed. Allowing it to be stack allocated, meaning that it is fast to do operations on it.
"""

# ╔═╡ c26c5949-e87a-4211-a49d-08a408ff3cdb
begin
	v1 = rand(3)
	v2 = rand(3)
	@benchmark v1 .+ v2
end

# ╔═╡ 109b282a-c156-48fe-a29d-ab543ec1a298
begin
	sv1 = SVector{3}(rand(3))
	sv2 = SVector{3}(rand(3))
	@benchmark sv1 .+ sv2
end

# ╔═╡ f489b8e0-1f5d-45f7-b5c3-03f09a5ef968
md"""
We also define our own class called a direction, which contains the zenith angle, `θ`, and azimuth angle, `ϕ`, as well as the projection onto the x-, y-, and z-axes. This last bit may not be necessary but I had it so I kept it. I might combine these into a SVector{3} which would probably be more efficient.
"""

# ╔═╡ c48ccde2-e02e-44dd-9b4a-400bae03b1ab
md"""
We are also able to define a `Box` which defines the region where we will consider injection. As long as the principle axes of a `Box` are aligned with our coordinate system, we can define a `Box` only by the two, opposite corners. Since this box will be the region in which our spline is defined, aligning the axes of the box and our coordinate system is a natural choice.
"""

# ╔═╡ 617a55f9-4d41-45a3-b386-04aeafef1f18
md"""
We may often want to know if points are inside or outside a `Box`
"""

# ╔═╡ 5e8a10d6-e5b6-4d48-9d46-45608b096d61
begin
	pt1 = SVector{3}([0, 10, 0])
	is_inside(pt1, b)
end

# ╔═╡ 14c3bea4-5ce7-4b14-a40b-830b247ac708
md"""
We make the actual Geometry, which has a spline of the Colca valley in the `valley` field, a `Box` defining the injection region in the the `box` field, and a `SVector{3}` defining where TAMBO is. This last attribute is mostly used for internal conversions and likely will not need to be interacted with.
"""

# ╔═╡ 59344208-d1b1-498f-ab00-56c67b2dd1c3
md"""
There are a handfull of ways to construct a `Geometry`. The first is only with a spline. In this case, we just put TAMBO in the middle of the spline, on surface of the mountain. This is a mostly just a convenience for developing.
"""

# ╔═╡ 13d9450f-a3ea-4849-b377-c440f54f9c0b
md"""
TAMBO is always at the center of the coordinate system, and since we placed it on the mountain, we expect the spline evaluated at the origin should be 0.
"""

# ╔═╡ 80705622-261a-4c5f-b8a6-adaa403c157c
geo1(100m, 100m) / m

# ╔═╡ c9e51ad8-f1dd-4b80-94b5-3797516513bd
md"""
Let's take a look at the mountain in this case !
"""

# ╔═╡ ee3fc122-b936-4c0e-8efc-affd42a4c84c
@bind azimuth1 Slider(0:90; show_value=true, default=15)

# ╔═╡ 9effba5c-bd68-4fb4-8368-e47df97c1035
@bind declination1 Slider(0:90 ; show_value=true, default=65)

# ╔═╡ 6464f79f-0e7c-4f2e-b2e1-8db1bfb87dc1
begin
	using Plots: palette, cgrad, surface, plot!, scatter3d!
	ax = surface(
    	LinRange(geo1.box.c1[1], geo1.box.c2[1], 100),
    	LinRange(geo1.box.c1[2], geo1.box.c2[2], 100),
    	geo1.valley,
    	c=cgrad(palette([
			:skyblue3,
			:skyblue2,
			:navajowhite3,
			:navajowhite3,
			:goldenrod4,
			:goldenrod4,
			:olivedrab,
			:olivedrab,
			:green,
			:green, 
			:green,
			:green
		])),
		colorbar=false,
	)
	plot!(ax, camera=(azimuth1, declination1))
	scatter3d!(ax, [0], [0], [0], markershape=:star, label="TAMBO", markersize=7)
end

# ╔═╡ 6cb06b77-33ff-48e6-a45d-abb1ff4d333a
md"""
Another way to construct a geometry is to specify only the xy-offset of TAMBO. In this case TAMBO will be placed on the mountain at this point.
"""

# ╔═╡ c15f9b51-2710-4a47-9eb9-5d064d2f87ac
begin
	xyoffset = SVector{2}([12500m, 12000m])
	geo2 = Geometry(spl, xyoffset)
end

# ╔═╡ be0950b1-9894-4db0-9b9f-12c77eae3c55
md"""
Once again, we placed the TAMBO offset on the mountain, and the offset is the origin of our coordinate system. Thus, we should expect that the valley evaluated at the center should be 0.
"""

# ╔═╡ 6ae70d11-e5cb-497a-80dc-134a539c3854
geo2(0,0)

# ╔═╡ ece09865-1b47-4c85-9f8c-53d49906d874
md"""
Now, let's take a peek at this valley.
"""

# ╔═╡ 17da74ea-3a82-41f5-8422-c88a52bf71b4
@bind azimuth2 Slider(0:90; show_value=true, default=0)

# ╔═╡ bedf04fa-fffd-468d-8e1b-23d90c0016e5
@bind declination2 Slider(0:90 ; show_value=true, default=65)

# ╔═╡ 2a457860-7221-4715-87ab-aacd8b84d68c
begin
	ax2 = surface(
    	LinRange(geo2.box.c1[1], geo2.box.c2[1], 100),
    	LinRange(geo2.box.c1[2], geo2.box.c2[2], 100),
    	geo2.valley,
    	c=cgrad(palette([
			:skyblue3,
			:skyblue2,
			:navajowhite3,
			:navajowhite3,
			:goldenrod4,
			:goldenrod4,
			:olivedrab,
			:olivedrab,
			:green,
			:green, 
			:green,
			:green
		])),
		colorbar=false,
	)
	plot!(ax2, camera=(azimuth2, declination2))
	scatter3d!(ax2, [0], [0], [0], markershape=:star, label="TAMBO", markersize=7)
end

# ╔═╡ 7ebf33b9-4aa9-4671-a920-d6886b8cb482
md"""
Lastly, we can tell the code directly where to put the TAMBO offset. This is useful in the case that the center of mass of the detector lies off the surface of the mountain because of some concavity or convexity.
"""

# ╔═╡ d5207b17-556f-43ed-b4f4-34c69b47cd3e
begin
	xyzoffset = SVector{3}([8000m, 12000m, 3000m])
	geo3 = Geometry(spl, xyzoffset)
end

# ╔═╡ 2808c118-d4d2-48fc-a162-a2fbae5bec4f
md"""
In this `Geometry` we have placed the offset suspended in the valley, and thus we would expect the valley, evaluated at (0,0), to have a negative value.
"""

# ╔═╡ e763ca9b-9118-4792-9362-44525b160605
geo3(0,0)/m

# ╔═╡ 6377ae93-1557-425f-9bad-46121d32779f
@bind azimuth3 Slider(0:90; show_value=true, default=0)

# ╔═╡ a356899f-fbe3-44f7-a981-e118bd3b6dc5
@bind declination3 Slider(0:90 ; show_value=true, default=65)

# ╔═╡ 112bcacd-4bb4-4147-9531-0136e3cb5f43
begin
	ax3 = surface(
    	LinRange(geo3.box.c1[1], geo3.box.c2[1], 100),
    	LinRange(geo3.box.c1[2], geo3.box.c2[2], 100),
    	geo3.valley,
    	c=cgrad(palette([
			:skyblue3,
			:skyblue2,
			:navajowhite3,
			:navajowhite3,
			:goldenrod4,
			:goldenrod4,
			:olivedrab,
			:olivedrab,
			:green,
			:green, 
			:green,
			:green
		])),
		colorbar=false,
	)
	plot!(ax3, camera=(azimuth3, declination3))
	scatter3d!(ax3, [0], [0], [0], markershape=:star, label="TAMBO", markersize=7)
end

# ╔═╡ Cell order:
# ╠═d260f8b0-d780-11ec-1e84-695ebd7b2393
# ╟─7ded89b6-496e-4bed-b3a2-53e8be42f6a3
# ╠═78fb3deb-4202-4e91-96b3-b89fbd99eaf4
# ╠═c26c5949-e87a-4211-a49d-08a408ff3cdb
# ╠═109b282a-c156-48fe-a29d-ab543ec1a298
# ╟─f489b8e0-1f5d-45f7-b5c3-03f09a5ef968
# ╠═4b006b1e-f0bc-4afc-8544-4b1d44d29afe
# ╟─c48ccde2-e02e-44dd-9b4a-400bae03b1ab
# ╠═0312f8e0-50b5-41d0-9342-1d0835c4f846
# ╟─617a55f9-4d41-45a3-b386-04aeafef1f18
# ╠═f28c2ca0-e241-4e02-bfc2-445688c81bc6
# ╠═5e8a10d6-e5b6-4d48-9d46-45608b096d61
# ╟─14c3bea4-5ce7-4b14-a40b-830b247ac708
# ╟─59344208-d1b1-498f-ab00-56c67b2dd1c3
# ╠═38ac33c5-a9e6-40cd-9244-afd6547f6488
# ╟─13d9450f-a3ea-4849-b377-c440f54f9c0b
# ╠═c04d426a-0e52-4d5b-a0a2-0b16c4688e93
# ╠═80705622-261a-4c5f-b8a6-adaa403c157c
# ╟─c9e51ad8-f1dd-4b80-94b5-3797516513bd
# ╠═9ea8bd84-efa1-4f54-a077-bf8faa9698d4
# ╟─ee3fc122-b936-4c0e-8efc-affd42a4c84c
# ╟─9effba5c-bd68-4fb4-8368-e47df97c1035
# ╟─6464f79f-0e7c-4f2e-b2e1-8db1bfb87dc1
# ╟─6cb06b77-33ff-48e6-a45d-abb1ff4d333a
# ╠═c15f9b51-2710-4a47-9eb9-5d064d2f87ac
# ╟─be0950b1-9894-4db0-9b9f-12c77eae3c55
# ╠═6ae70d11-e5cb-497a-80dc-134a539c3854
# ╟─ece09865-1b47-4c85-9f8c-53d49906d874
# ╟─17da74ea-3a82-41f5-8422-c88a52bf71b4
# ╟─bedf04fa-fffd-468d-8e1b-23d90c0016e5
# ╠═2a457860-7221-4715-87ab-aacd8b84d68c
# ╟─7ebf33b9-4aa9-4671-a920-d6886b8cb482
# ╠═d5207b17-556f-43ed-b4f4-34c69b47cd3e
# ╟─2808c118-d4d2-48fc-a162-a2fbae5bec4f
# ╠═e763ca9b-9118-4792-9362-44525b160605
# ╟─6377ae93-1557-425f-9bad-46121d32779f
# ╟─a356899f-fbe3-44f7-a981-e118bd3b6dc5
# ╟─112bcacd-4bb4-4147-9531-0136e3cb5f43
