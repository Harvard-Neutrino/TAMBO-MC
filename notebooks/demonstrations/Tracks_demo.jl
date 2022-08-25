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

# ╔═╡ 2eab217c-e0ff-11ec-0764-c16da7ed9296
begin
	using Pkg
	using StaticArrays
	Pkg.activate("/Users/jlazar/research/TAMBO-MC/tambo/")
	tambo_path = abspath("../../tambo/src/")
	using Tambo
end

# ╔═╡ 17cfa11b-041c-44ab-8361-ed99fe068a03
using PlutoUI: Slider

# ╔═╡ eba4251d-d763-4bd7-9964-a221c476fdfb
md"""
Next, we move on to `Tracks.jl`. Note that `Geometries.jl` is included within `Tracks.jl`. Thus, even though we will be using some `struct`s from it, we do not need to `include` it again.
"""

# ╔═╡ 8091543e-e591-460b-9daf-ac3ca59440ed
md"""
This file is mainly the `Track` `struct` with some functions for calculating qunatities of interest for the this. As you can see there are four different ways to construct a `Track`. We will give an example of each in the following.
"""

# ╔═╡ 96d8aada-9ef1-40c8-9a04-0484026d67c5
methods(Track)

# ╔═╡ 4887d433-b57f-4b09-9f49-66ddabcebcda
md"""
First, if we pass in a single point, this will be taken as the end point of the track, and the starting point will be the origin.
"""

# ╔═╡ dccb675d-0fd6-409b-a528-4e6ce51991b5
begin
	end_point = SVector{3}([1, 2, 3])
	t1 = Tambo.Track(end_point)
end

# ╔═╡ 99f23374-b1bc-4597-9b4c-4964f981236d
md"""
Next, if we pass in two points, the first point will be taken to be the starting point, and the second will be the end point.
"""

# ╔═╡ 9b68b3ea-f68a-4aaf-b0b4-5b40ee268460
begin
	init_point = SVector{3}([-1.0, -2.0, -3.0])
	t2 = Tambo.Track(init_point, end_point)
end

# ╔═╡ a9c492ba-7611-48db-b9a5-af2be86b07ca
md"""
Lastly, we can define a track by specificying an initial point, initial direction, and requiring that the track end on a bounding `Box`. This will be of particular importance to us since the bounds of the TAMBO simulation are a `Box`.
"""

# ╔═╡ 519d9549-a817-491a-b37c-851b8699f266
begin
	d = Tambo.Direction(π/2, 13π/12)
	spl = load_spline(tambo_path * "../../resources/tambo_spline.jld2")
	geo = Tambo.Geometry(spl)
	t3 = Tambo.Track(init_point, d, geo.box)
end

# ╔═╡ a74ae89e-250d-4620-b17c-45c86d329b12
md"""
Let's make sure that the end point of the `Track` is in fact on one of the `Box` boundaries.
"""

# ╔═╡ e88642d9-8e2a-4a7c-ab09-a608ad111ef3
any(geo.box.c1 .== t3.fpoint) || any(geo.box.c2 .== t3.fpoint)

# ╔═╡ 192eaa70-5767-4133-b60e-efdaedd8c1c3
md"""
We can find the position of the track at any point along the trajectory by calling the `Track` as if it were a function.
"""

# ╔═╡ b4fe113f-fafa-450d-bb5e-20c26cf7a616
t3(0.5) ./ units[:m]

# ╔═╡ b67bf248-2329-4ffd-8673-25ef19b09029
md"""
The argument of this function is the affine parameter that tells us how far along the track we are. It ALWAYS runs between 0 and 1, with 0 corresonding to the start of the `Track`, and 1 the end of the `Track`.
"""

# ╔═╡ d23f4f66-2d5d-4ca1-b3a9-d87c9d59719d
all(t3(0.0) .== t3.ipoint) && all(t3(1.0) .== t3.fpoint)

# ╔═╡ cf8c2add-9c72-4483-980e-ff512c0b160d
md"""
We can the physical distance the track has gone by multiplying the affine parameter by the norm of the `Track`
"""

# ╔═╡ 1c7496e2-2ca0-4576-b1f9-f987fda843b3
begin
	λ0 = 0.6
	λ0 * t3.norm == sqrt(sum((t3(λ0) .- t3.ipoint).^2))
end

# ╔═╡ 19446848-4195-4884-b388-0c30ea8d39d7
md"""
We can also reverse a `Track`, swapping the `ipoint` and `fpoint`. Let's confirm that the directions are (approximately) reversed.
"""

# ╔═╡ 1f3979b0-e281-45d3-8b98-664b2c38fe63
begin 
	t3r = reverse(t3)
	all(t3r.direction.proj .≈ -t3.direction.proj)
end

# ╔═╡ f5851bda-141b-4137-8060-31785d4f00bc
md"""
Now, let's take a look at what happens when we plop a `Track` into the TAMBO `Geometry`. First, we can calculate the total column depth.
"""

# ╔═╡ 16e7e5e7-a054-4aaf-95f5-86fb707318db
begin
	Tambo.totalcolumndepth(t3, geo) / units[:mwe]
end

# ╔═╡ 5b81ea74-62af-4b63-8dc9-f94195492e13
md"""
We can also watch it fly through the valley.
"""

# ╔═╡ 559860f2-1bf7-4553-aa34-a7f868a6e7ac
md"""
Viewing angles
"""

# ╔═╡ aa9eb026-0ce4-4c21-835a-674afa42ad88
@bind azimuth1 Slider(0:90; show_value=true, default=15)

# ╔═╡ 82c01d9b-321d-4232-855f-3cec3cc3f6a9
@bind declination1 Slider(0:90 ; show_value=true, default=65)

# ╔═╡ 9df7eca0-9268-4e64-b32f-e38e298b3e09
md"""
`Track` direction
"""

# ╔═╡ 7943180c-b01a-4ad4-a78d-c0d0c09a5847
@bind zenith2 Slider(0:180, show_value=true, default=85)

# ╔═╡ 339acc0b-5288-46e3-9ae3-cf5bf2aaba9f
@bind azimuth2 Slider(0:360, show_value=true, default=175)

# ╔═╡ 9cefadb8-8fd3-4d90-9596-a1b06838ec3b
begin
	d1 = Tambo.Direction(zenith2*pi/180, azimuth2*pi/180)
	t4 = Tambo.Track(init_point, d1, geo.box)
	md"Affine parameter"
end

# ╔═╡ ab0e8923-84be-4f8b-ac20-e84d710f6de8
@bind λcut Slider(LinRange(0, 1, 101), default=1.0, show_value=true)

# ╔═╡ 08448d04-3ebf-4fdc-ba17-76acaabea827
begin
	using Plots: palette, cgrad, surface, plot!, scatter3d!, plot3d!
	ax = surface(
    	LinRange(geo.box.c1[1], geo.box.c2[1], 100),
    	LinRange(geo.box.c1[2], geo.box.c2[2], 100),
    	geo.valley,
    	c=cgrad(palette([
        :skyblue3,
        :skyblue2,
        :skyblue2,
        :skyblue2,
        :green, 
        :chartreuse4,
        :chartreuse4,
        :goldenrod4,
        :goldenrod4,
        :goldenrod4,
        :goldenrod4,
        :navajowhite3,
        :white,
        :white,
        :white,
        :white,
        :white,
	    ])),
		colorbar=false,
	)
	plot!(ax, camera=(azimuth1, declination1), zlim=(-1e10, 2.53387e10))
	λλ = LinRange(0, 1, 101)
	λλ = λλ[λλ.<=λcut]
	t = t4.(λλ)
	plot3d!(
		ax,
		getindex.(t ,1 ), 
		getindex.(t, 2 ),
		getindex.(t, 3 ),
		line_z=λλ,
		c=cgrad(:rainbow),
		clim=(0,1),
		lw=2,
		legend=false
	)
	
end

# ╔═╡ ce8390d6-83b4-4fef-965a-6bfc7029c40e
md"""
Is the `Track` inside the mountain ?
"""

# ╔═╡ 0ac14b31-28db-4829-9441-657f7fb42212
Tambo.inside(t4(λcut), geo.valley)

# ╔═╡ a7a120e9-d2a3-4d78-95a0-5841a14bf884
md"""
Column depth so far (mwe)
"""

# ╔═╡ a9aaf6ce-2759-4611-8018-9e7eaf2a882c
begin
	Tambo.columndepth(t4, λcut, geo) / units[:mwe]
end

# ╔═╡ 39419e4e-9812-4cbb-bbf3-d4e0664993b4


# ╔═╡ 734c4fc7-98cf-4344-8b33-3f27d705c617


# ╔═╡ 950def7e-1a76-470c-8609-64ae95ccaea7


# ╔═╡ 04580e99-0cd4-42f3-9806-1e049cf0c15e


# ╔═╡ 2365bc95-f622-427b-8ad1-537de52827fd


# ╔═╡ Cell order:
# ╟─eba4251d-d763-4bd7-9964-a221c476fdfb
# ╠═2eab217c-e0ff-11ec-0764-c16da7ed9296
# ╟─8091543e-e591-460b-9daf-ac3ca59440ed
# ╟─96d8aada-9ef1-40c8-9a04-0484026d67c5
# ╟─4887d433-b57f-4b09-9f49-66ddabcebcda
# ╠═dccb675d-0fd6-409b-a528-4e6ce51991b5
# ╟─99f23374-b1bc-4597-9b4c-4964f981236d
# ╠═9b68b3ea-f68a-4aaf-b0b4-5b40ee268460
# ╟─a9c492ba-7611-48db-b9a5-af2be86b07ca
# ╠═519d9549-a817-491a-b37c-851b8699f266
# ╟─a74ae89e-250d-4620-b17c-45c86d329b12
# ╠═e88642d9-8e2a-4a7c-ab09-a608ad111ef3
# ╟─192eaa70-5767-4133-b60e-efdaedd8c1c3
# ╠═b4fe113f-fafa-450d-bb5e-20c26cf7a616
# ╟─b67bf248-2329-4ffd-8673-25ef19b09029
# ╠═d23f4f66-2d5d-4ca1-b3a9-d87c9d59719d
# ╟─cf8c2add-9c72-4483-980e-ff512c0b160d
# ╠═1c7496e2-2ca0-4576-b1f9-f987fda843b3
# ╟─19446848-4195-4884-b388-0c30ea8d39d7
# ╠═1f3979b0-e281-45d3-8b98-664b2c38fe63
# ╟─f5851bda-141b-4137-8060-31785d4f00bc
# ╠═16e7e5e7-a054-4aaf-95f5-86fb707318db
# ╟─5b81ea74-62af-4b63-8dc9-f94195492e13
# ╠═17cfa11b-041c-44ab-8361-ed99fe068a03
# ╟─559860f2-1bf7-4553-aa34-a7f868a6e7ac
# ╟─aa9eb026-0ce4-4c21-835a-674afa42ad88
# ╟─82c01d9b-321d-4232-855f-3cec3cc3f6a9
# ╟─9df7eca0-9268-4e64-b32f-e38e298b3e09
# ╟─7943180c-b01a-4ad4-a78d-c0d0c09a5847
# ╟─339acc0b-5288-46e3-9ae3-cf5bf2aaba9f
# ╟─9cefadb8-8fd3-4d90-9596-a1b06838ec3b
# ╟─ab0e8923-84be-4f8b-ac20-e84d710f6de8
# ╟─ce8390d6-83b4-4fef-965a-6bfc7029c40e
# ╟─0ac14b31-28db-4829-9441-657f7fb42212
# ╟─a7a120e9-d2a3-4d78-95a0-5841a14bf884
# ╟─a9aaf6ce-2759-4611-8018-9e7eaf2a882c
# ╠═08448d04-3ebf-4fdc-ba17-76acaabea827
# ╠═39419e4e-9812-4cbb-bbf3-d4e0664993b4
# ╠═734c4fc7-98cf-4344-8b33-3f27d705c617
# ╟─950def7e-1a76-470c-8609-64ae95ccaea7
# ╠═04580e99-0cd4-42f3-9806-1e049cf0c15e
# ╠═2365bc95-f622-427b-8ad1-537de52827fd
