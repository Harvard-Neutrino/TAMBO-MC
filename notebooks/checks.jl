### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 505936e8-592b-11ed-0134-f90a549c41ff
begin
	using PlutoUI: Slider
	using Plots
	tambo_path = abspath("../Tambo/")
	using Pkg
	Pkg.activate(tambo_path)
	using Tambo
	using StaticArrays
	using JLD2
	units = Tambo.units
end;

# ╔═╡ 7640a98f-eb75-4101-8bc4-ac91e19becc5
Pkg.activate()


# ╔═╡ 64056cbe-280a-4b38-8ec4-711c3d0f9233
Pkg.add("WebIO")

# ╔═╡ 110bf951-1466-4152-aa99-22d7d7ed94da
using PlotlyJS

# ╔═╡ 77b08741-1287-4802-969b-abde2e78cf69
begin
		spl = Tambo.load_spline(joinpath(tambo_path, "../resources/tambo_spline.jld2"))

	latlongmin = Tambo.mincoord_fromfile(
		filename=realpath("../resources/ColcaValleyData.txt")
	)
	mine_x, mine_y = Tambo.latlong_to_xy(Tambo.minesite_coord, latlongmin)
	test_x, test_y = Tambo.latlong_to_xy(Tambo.testsite_coord, latlongmin)
	tapay_x = 42.48869264279284 * units.km - mine_x
	tapay_y = 17.927740080846274 * units.km - mine_y
	xyoffset = SVector{2}([mine_x, mine_y])
	geo = Tambo.Geometry(spl, xyoffset)
end

# ╔═╡ abe7f8a2-0b12-491c-9995-8fc486326cbe
@bind azimuth Slider(0:90; show_value=true, default=15)

# ╔═╡ d67f5895-720f-47b0-a826-35450e90111b
@bind declination Slider(0:90; show_value=true, default=15)

# ╔═╡ 10a1e9d9-fe9d-43ea-a86d-26bd0ca0b6eb
begin
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
	n2 = 1000
	x2 = LinRange(geo.box.c1[1], geo.box.c2[1], n2)
	y2 = LinRange(geo.box.c1[2], geo.box.c2[2], n2)
	ax2 = surface(
    	x2 ./ units.km,
    	y2./ units.km,
    	reshape([geo.valley(x,y)/units.km for x in x2 for y in y2], (n2,n2)),
    	c=cmap,
		colorbar=false,
	)
	plot!(ax2, camera=(azimuth, declination))
	scatter3d!(ax2, [0], [0], [0], markershape=:star, label="Mine site", markersize=7)
	scatter3d!(ax2, [(test_x - mine_x) / units.km], [(test_y - mine_y) / units.km], [0], markershape=:star, label="Test site", markersize=7, color=:gold)
	scatter3d!(ax2, [tapay_x/units.km], [tapay_y/units.km], [0], label="Tapay",
		color="blue")
end

# ╔═╡ 894f2d13-e0e1-4d06-b680-257e06eb4fe3
begin
	injector = Simulator()
end

# ╔═╡ df4419e0-dd32-409a-9742-ad40a6b9cd6c
f = jldopen("../injected_events")

# ╔═╡ ccb41fce-3b22-4619-8b89-0236033799c5
events = f["injected_events"]

# ╔═╡ f9e96d4c-b5e0-4692-8452-cc9e56eacf8a
@bind azimuth1 Slider(0:90; show_value=true, default=15)

# ╔═╡ a2c83097-c8ea-40d1-933c-66981de9204e
@bind declination1 Slider(0:90; show_value=true, default=15)

# ╔═╡ 7a7c40c9-70d0-4fb7-a713-e0c082a9e780
begin
	xs = LinRange(-20, 20, 100) * units.km
	ys = LinRange(-16, 16, 100) * units.km
	sct = PlotlyJS.scatter(
		x=getindex.(events["final_state"]["position"][1:10_000], 1) ./ units.km,
		y=getindex.(events["final_state"]["position"][1:10_000], 2) ./ units.km,
		# getindex.(events["final_state"]["position"], 3) ./ units.km,
		# alpha=0.1,
		# label=""
	)
	# plot!(
	# 	sct,
	# 	contour(
	# 		x=xs / units.km, y=ys / units.km,
	# 		z=reshape([geo.valley(x, y) for x in xs for y in ys], (100,100))
	# 	)
	# )
	# display(sct)
end

# ╔═╡ 3e9dcfda-1992-49e8-a3d0-c5765dc45c0e
PlotlyJS.plot(PlotlyJS.contour(
    x=[-9, -6, -5 , -3, -1], # horizontal axis
    y=[0, 1, 4, 5, 7], # vertical axis
    z=[
        10      10.625      12.5       15.625     20
        5.625    6.25       8.125      11.25      15.625
        2.5      3.125      5.         8.125      12.5
        0.625    1.25       3.125      6.25       10.625
        0        0.625      2.5        5.625      10
    ]'
))

# ╔═╡ 349448d4-eedb-4b77-8f21-6d2311103c33
methods(PlotlyJS.scatter)

# ╔═╡ Cell order:
# ╠═505936e8-592b-11ed-0134-f90a549c41ff
# ╠═77b08741-1287-4802-969b-abde2e78cf69
# ╠═abe7f8a2-0b12-491c-9995-8fc486326cbe
# ╠═d67f5895-720f-47b0-a826-35450e90111b
# ╠═10a1e9d9-fe9d-43ea-a86d-26bd0ca0b6eb
# ╠═894f2d13-e0e1-4d06-b680-257e06eb4fe3
# ╠═df4419e0-dd32-409a-9742-ad40a6b9cd6c
# ╠═ccb41fce-3b22-4619-8b89-0236033799c5
# ╟─f9e96d4c-b5e0-4692-8452-cc9e56eacf8a
# ╟─a2c83097-c8ea-40d1-933c-66981de9204e
# ╠═7a7c40c9-70d0-4fb7-a713-e0c082a9e780
# ╠═110bf951-1466-4152-aa99-22d7d7ed94da
# ╠═7640a98f-eb75-4101-8bc4-ac91e19becc5
# ╠═64056cbe-280a-4b38-8ec4-711c3d0f9233
# ╠═3e9dcfda-1992-49e8-a3d0-c5765dc45c0e
# ╠═349448d4-eedb-4b77-8f21-6d2311103c33
