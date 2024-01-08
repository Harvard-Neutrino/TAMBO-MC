### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ 98e9ed62-4f03-11ed-2672-67d543b60dfe
begin
	tambo_path = abspath("../../Tambo/")
	using Pkg
	Pkg.activate(tambo_path)
	using Tambo
	units = Tambo.units;
end

# ╔═╡ 54c35ae0-070d-448c-929b-ef194528972a
begin
	spl = Tambo.load_spline(joinpath(tambo_path, "../resources/tambo_spline.jld2"))

	using StaticArrays: SVector
	using Plots: palette, cgrad, surface, plot!, scatter3d!
	cmap = cgrad(palette([
			 :skyblue3,
			:skyblue2,
			:skyblue2,
			:skyblue2,
			:chartreuse4,
			:chartreuse4,
			:chartreuse4,
			:goldenrod4,
			:goldenrod4,
			:goldenrod4,
			:peru,
			:peru,
			:peru,
			:navajowhite3,
		]))
	# proposed_x = 39.930033123220795 * units.km
	# proposed_y = 17.090219893359233 * units.km
	# proposed_x2 = 20.65344955402193 * units.km- proposed_x
	# proposed_y2 = 9.159464140289597 * units.km - proposed_y
	# tapay_x = 42.48869264279284 * units.km - proposed_x
	# tapay_y = 17.927740080846274 * units.km - proposed_y
	# xyoffset = SVector{2}([proposed_x, proposed_y])
	geo = Tambo.Geometry(spl)
	n = 200
	x = LinRange(geo.box.c1[1], geo.box.c2[1], n)
	y = LinRange(geo.box.c1[2], geo.box.c2[2], n)
	ax = surface(
    	x ./ units.km,
    	y./ units.km,
    	reshape([geo.valley(x,y)/units.km for x in x for y in y], (n, n)),
    	c=cmap,
		colorbar=false,
	)
	plot!(ax, camera=(0, 0))
	scatter3d!(ax, [0], [-13], [0], markershape=:star, label="Proposed site 1", markersize=7)
	# scatter3d!(ax, [proposed_x2 / units.km], [proposed_y2 / units.km], [0], markershape=:star, label="Proposed site 2", markersize=7, color=:gold)
	# scatter3d!(ax, [tapay_x/units.km], [tapay_y/units.km], [0], label="Tapay", color="blue")
end

# ╔═╡ 7196e85f-4d27-455e-af63-97305f3070b1
begin
	using Plots
	kwargs = Dict(
		# :xlim=>(minimum((geo.box.c1[1], geo.box.c2[1])), maximum((geo.box.c1[1], geo.box.c2[1]))),
	 #    :ylim=>(minimum((geo.box.c1[2], geo.box.c2[2])), maximum((geo.box.c1[2], box.c2[2]))),
	 #    :zlim=>(minimum((geo.box.c1[2], geo.box.c2[3])), maximum((geo.box.c1[3], geo.box.c2[3]))),
	    :st=>:surface,
	    :alpha=>0.9,
	    :c=>cmap,
	    :colorbar=>false,
	    :legend=>false,
	)
	xx = LinRange(geo.box.c1[1], geo.box.c2[1], 100)
	yy = LinRange(geo.box.c1[2], geo.box.c2[2], 100)
	l = Plots.@layout [a{0.5w}  c ; d]
	plt = surface(
    	x ./ units.km,
    	y./ units.km,
    	reshape([geo.valley(x,y)/units.km for x in xx for y in yy], (100, 100)),
    	c=cmap,
		colorbar=false,
		layout=l,
		
		camera=(0, 90);
		kwargs...
	)
	# This should be done with an imshow-like thing.... I think
	surface!(
		plt[2],
		x ./ units.km,
    	y./ units.km,
    	reshape([geo.valley(x,y)/units.km for x in x for y in y], (n, n)),
    	c=cmap,
		colorbar=false,
		camera=(0, 90);
		kwargs...
	)
	
	# if show_track
	#     λλ  = LinRange(0, 1, 50)
	#     pts = [tr.ipoint+λ*tr.direction for λ in λλ]
	#     segments = split_for_plot([pts[idx] for idx in 1:3]...)
	#     for segment in segments
	#         for idx in 1:3
	#             plot!(plt[idx], segment[1:3]..., color=segment[4], legend=false)
	#         end
	#     end
	# end
	
	# display(plt)
	# Plots.savefig(plt, "/Users/jlazar/research/TAMBO-MC/figures/fake_valley.pdf")
	
	# if gifify && ~show_track # These options don't play well together. Sorry
	# trayectory = Trayectory(tr, dλ)    
	#     @gif while trayectory.current_λ <=1
	#         step!(trayectory)
	#         segments = split_for_plot(trayectory)
	#         for segment in segments
	#             for idx in 1:3
	#                 plot!(plt[idx], segment[1:3]..., color=segment[4], legend=false)
	#             end
	#         end
	#         plot!(plt[1], camera=(minimum((90, trayectory.current_λ*100)), minimum((90, trayectory.current_λ*100))))
	#     end every 1
	# end
end

# ╔═╡ f384e6a1-69ce-4d1d-b81b-bf93c3ec5fbb
begin
	mu_injector = Injector()
	mu_injector.ν_pdg = 13
	tau_injector = Injector()
	injected_mus = mu_injector()
	injected_taus = tau_injector()
end

# ╔═╡ 5d3c38ba-d8da-445c-8885-7e52ce35c6e3
begin
	p = SVector{3}([0, -13*Tambo.units.km, 0])
	d = Tambo.Direction(1.0, -1.0, 0.0)
	mu = Tambo.Particle(13, injected_mus[1].final_state.energy, p, d, nothing)
	tau = Tambo.Particle(15, injected_mus[1].final_state.energy, p, d, nothing)
	propped_mu = Tambo.propagate(mu, geo)
	propped_tau = Tambo.propagate(tau, geo)
end

# ╔═╡ 78806267-fe52-44ab-93bd-5293b8dfea47
begin
	t1 = Tambo.Track(p, d, geo.box)
	t2 = Tambo.Track(p, reverse(d), geo.box)
	t = Tambo.Track(t1.fpoint, t2.fpoint)
	f(λ) = t(λ).x - p.x
	using Roots: find_zero
	find_zero(f, 0.5)
end

# ╔═╡ 6cd6aadb-9d58-4308-b7d6-dc15f91dfa39
t1

# ╔═╡ 3ae0189a-e409-4235-968d-6c77db720757
l

# ╔═╡ Cell order:
# ╠═98e9ed62-4f03-11ed-2672-67d543b60dfe
# ╠═6cd6aadb-9d58-4308-b7d6-dc15f91dfa39
# ╠═54c35ae0-070d-448c-929b-ef194528972a
# ╟─f384e6a1-69ce-4d1d-b81b-bf93c3ec5fbb
# ╠═5d3c38ba-d8da-445c-8885-7e52ce35c6e3
# ╠═78806267-fe52-44ab-93bd-5293b8dfea47
# ╠═3ae0189a-e409-4235-968d-6c77db720757
# ╠═7196e85f-4d27-455e-af63-97305f3070b1
