### A Pluto.jl notebook ###
# v0.19.26

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

# ╔═╡ 75cc79a1-a15e-4631-a275-3ee3053aaa5d
begin 
	using Dierckx
	tambopath = realpath("/Users/pavelzhelnin/Documents/physics/TAMBO/tambo");
	using Pkg; 
	Pkg.activate(tambopath)
	
	using Tambo
	simulator = SimulationConfig("/Users/pavelzhelnin/Documents/physics/TAMBO/notebooks/WhitePaper_EarthMantle_test.jld2")
	geo = Geometry("$(tambopath)/../resources/tambo_spline.jld2", simulator.tambo_coordinates);
end 

# ╔═╡ c8b8d3bd-dedf-4333-b0e7-764e7aa2e9cd
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

# ╔═╡ 3ac8cead-b089-4720-91cc-10fef3dbcd5a
@bind λcut Slider(LinRange(0, 1, 101), default=1.0, show_value=true)

# ╔═╡ 162daaf3-14ba-4099-9308-8f2c7257eede


# ╔═╡ 60518038-b4d8-4000-bd9e-d77015217cab
@bind azimuth1 Slider(-90:90; show_value=true, default=15)

# ╔═╡ b08b4b48-4e39-4a3b-858c-47a2666c3b04
function full_particle_viewer(index)
	
		decay_state = simulator["proposal_events"]["decay_products"][index][1]
		state = simulator["proposal_events"]["propped_state"][index]
		injected_event = simulator["injected_events"][index]
		CC_pos = injected_event.final_state.position/units.m
		taudecay_pos = decay_state.position/units.m
		zenith = state.direction.θ*180/pi
    	azimuth = state.direction.ϕ*180/pi

	df = CSV.File("/Users/pavelzhelnin/Documents/physics/TAMBO/notebooks/WhitePaper_EarthMantle.csv")

	result = index in df["index"]


		track = Tambo.Track(CC_pos,taudecay_pos)

		xx = LinRange(400 + geo.box.c1[1]/units.m, geo.box.c2[1]/units.m, 100)
		yy = LinRange(400 + geo.box.c1[2]/units.m, geo.box.c2[2]/units.m, 100)

		kwargs = Dict(
          :st=>:surface,
          :alpha=>0.9,
          :c=>cgrad(palette([:skyblue3, :skyblue2, :navajowhite3, :navajowhite3, :goldenrod4, :goldenrod4, :olivedrab, :olivedrab, :green, :green, :green, :green])),
          :colorbar=>false,
          :legend=>false,
         )
	f(x,y) = geo(x*units.m,y*units.m)/units.m

	l   = @layout [a{0.9w} ; b{0.5w} c{0.5w}]
	plt = surface(xx, yy, f, layout = l, axis=([], false), camera = (azimuth1,0); kwargs...)

	surface!(plt[3],xx, yy, f, layout = l, axis=([], false), camera = (azimuth1,90); kwargs...)


	plane = Tambo.Plane(minesite_normal_vec,minesite_coord,geo)

	p_phi = plane.n̂.ϕ
    p_theta = plane.n̂.θ 

    zdir = cos(p_theta)
    ydir = sin(p_theta) * sin(p_phi)
    xdir = sin(p_theta) * cos(p_phi)

	plane_x = LinRange(-5000, 5000, 100)

	pp(x,z) = (((ydir)*x)+(z*zdir))/(xdir)
	
	print("this is tau decay position")
	println(taudecay_pos)
	print("this is the CC interaction position")
	println(CC_pos)
	scatter!(plt[1],[taudecay_pos.x],[taudecay_pos.y],[taudecay_pos.z],color="red")
    scatter!(plt[1],[CC_pos.x],[CC_pos.y],[CC_pos.z],color="green")
	scatter!(plt[1],[0],[0],[0],color="blue")

	scatter!(plt[3],[taudecay_pos.x],[taudecay_pos.y],[taudecay_pos.z],color="red")
    scatter!(plt[3],[CC_pos.x],[CC_pos.y],[CC_pos.z],color="green")
	scatter!(plt[3],[0],[0],[0],color="blue")

	if result == true

		ind = findfirst(x -> x == index, df["index"])
		x_inter = (df["inter_x"][ind].*1000)
    	y_inter = (df["inter_y"][ind].*1000)
    	z_inter = (df["inter_z"][ind].*1000)
		scatter!(plt[1],[x_inter],[y_inter],[z_inter],color="black")

		
		track = Tambo.Track(CC_pos,SVector{3,Float64}(x_inter, y_inter, z_inter))
		println("This is the point of interception on TAMBO plane: $x_inter, $y_inter, $z_inter")
	else
		track = Tambo.Track(CC_pos,taudecay_pos)
	end 

	

	λλ = LinRange(0, 1, 101)
	λλ = λλ[λλ.<=λcut]
	t = track.(λλ)
	
	plot3d!(plt[1],
		getindex.(t ,1 ), 
		getindex.(t, 2 ),
		getindex.(t, 3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=false,
		legend=false
	)

	plot3d!(plt[3],
		getindex.(t ,1 ), 
		getindex.(t, 2 ),
		getindex.(t, 3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=false,
		legend=false
	)
	
	plot!(
		plt[2],
		λλ, 
		getindex.(t,3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=false,
		legend=false
	)
	
	bool1 = Tambo.inside(taudecay_pos,f)
	bool2 = Tambo.inside(CC_pos,f)
	println("Does the tau decay inside of the mountain: $bool1")
	println("Does the CC interaction happen inside of the mountain: $bool2")


	g_λλ = LinRange(-0.75, 1.75, 101)
	geo_track = track.(g_λλ)
	
	z_geo = [geo(ts[1]*units.m,ts[2]*units.m)/units.m for ts in geo_track]

	plot!(plt[2],g_λλ,z_geo,color = "green")
	scatter!(plt[2],[0],[getindex.(track.(0),3)],color="green")

	if result == true 
		scatter!(plt[2],[1],[getindex.(track.(1),3)],color="black")

		ind = findfirst(x -> x > taudecay_pos.z, getindex.(t,3))
		scatter!(plt[2],[λλ[ind]],[getindex.(t,3)[ind]],color="red")
	else 
		scatter!(plt[2],[1],[getindex.(track.(1),3)],color="red")
	end

	

	#plot!(plt[3])
	#scatter!(plt[3],[0],[getindex.(track.(0),3)],color="red")
	#scatter!(plt[3],[1],[getindex.(track.(1),3)],color="green")
	
end

# ╔═╡ bba917d0-f929-4ed7-ad5d-4f08912d7286
@bind declination1 Slider(0:90 ; show_value=true, default=65)

# ╔═╡ 2889faad-684d-4e93-9fc4-2d8cb07fc7a2
function particle_viewer_1d(index)
	
	decay_state = simulator["proposal_events"]["decay_products"][index][1]
	state = simulator["proposal_events"]["propped_state"][index]
	injected_event = simulator["injected_events"][index]
	CC_pos = injected_event.final_state.position/units.m
	taudecay_pos = decay_state.position/units.m
	zenith = state.direction.θ*180/pi
    azimuth = state.direction.ϕ*180/pi
	track = Tambo.Track(CC_pos,taudecay_pos)

	λλ = LinRange(0, 1, 101)
	λλ = λλ[λλ.<=λcut]
	t = track.(λλ)

	f(x,y) = geo(x*units.m,y*units.m)/units.m
	
	plot(
		λλ, 
		getindex.(t,3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=true,
		legend=false
	)
	bool1 = Tambo.inside(taudecay_pos,f)
	bool2 = Tambo.inside(CC_pos,f)
	println("Does the tau decay inside of the mountain: $bool1")
	println("Does the CC interaction happen inside of the mountain: $bool2")

	
	z_geo = [geo(ts[1]*units.m,ts[2]*units.m)/units.m for ts in t]

	plot!(λλ,z_geo,color = "green")
	scatter!([0],[getindex.(track.(0),3)],color="red")
	scatter!([1],[getindex.(track.(1),3)],color="green")

end

# ╔═╡ 0a203549-8d17-43c1-91c5-1e6018544d3f
particle_viewer_1d(11010)

# ╔═╡ 220345e2-83b3-4434-bcb4-3b5f5352b8bf
begin 
	state = simulator["proposal_events"]["propped_state"]

	zenith = [state.direction.θ*180/pi for state in state]

	histogram(zenith,bins=20,xlabel="Zenith angle in degrees")

end

# ╔═╡ b59160aa-5245-4827-ab0b-c6f73cc05402
begin 
	#state = simulator["proposal_events"]["propped_state"]

	azimuth = [state.direction.ϕ*180/pi for state in state]

	histogram(azimuth,bins=20,xlabel="Azimuthal angle in degrees")

end

# ╔═╡ e64ab6a5-edb6-482d-ab70-0aeeeb425431
begin 
	cos_zenith = [cos(state.direction.θ) for state in state]

	histogram(cos_zenith,bins=20,xlabel="cos(θ)")

end

# ╔═╡ 8c585aae-263a-4ba8-a8bb-4a3afb5f00be
begin 
decay_state = simulator["proposal_events"]["decay_products"]

taudecay_z = []
for ds in decay_state
	push!(taudecay_z,ds[1].position.z/units.m)
end
println(minimum(taudecay_z))
histogram([taudecay_z], bins=1000 ,xlabel="tau decay elevation (TAMBO) (m)",xlim=[-10000,12000],legend=false)
end 

# ╔═╡ e57cd392-fc79-4359-b461-d795b82ec082
begin 
propped_state = simulator["proposal_events"]["propped_state"]

taudecay_pos = []

for ds in decay_state
	push!(taudecay_pos,ds[1].position/units.m)
end


CC_pos = [state.position/units.m for state in propped_state]
	
# Initialize an array to store the distances
distances = []

# Compute distances between corresponding points in the sets
for (point1, point2) in zip(CC_pos, taudecay_pos)
    distance = norm(point1 - point2)
    push!(distances, distance)
end

energies = [state.energy/units.GeV for state in propped_state]

println(maximum(energies))
println(minimum(energies))

bar(energies, distances, bins=10 ,ylabel="L(E) meters", xlabel="energies")#xlim=(9e5,5e6),ylim=(0,200))
end 

# ╔═╡ 25f16065-881f-4d85-9091-70813bb6811b
bins = range(minimum(energies),maximum(energies),100)

# ╔═╡ b79157b5-9f07-413c-96cc-cd2905dfd86d
print(findfirst(x -> x < 500, bins))

# ╔═╡ fe9c6cea-463b-4474-bcc2-ddeb30b19391
inds = [findfirst(x->x > energy,bins) for energy in energies]

# ╔═╡ 58a0ddb1-d1b1-4c87-8e06-5d106df268f4
new_inds = [findall(x->x == ind, inds) for ind in range(1,100)]

# ╔═╡ c059d526-77dd-49c7-ab2b-7aa8cf9dbe69
begin
	median_lengths = []
	for ind in new_inds[1:100]
		try
		push!(median_lengths,median(distances[ind]))
		catch 
		println("No data")
		push!(median_lengths,0)
		end
	end
end

# ╔═╡ 9d746940-81ba-4b66-b4f0-e3b843a60cdd
bar(bins,median_lengths,xlabel="Energy (GeV)", ylabel = "L(E) (m)",legend=false)

# ╔═╡ 1df771f5-50e7-410b-aa85-a5c138c3d7a4


# ╔═╡ ab37cf29-61fe-4f8d-a251-5453e3e4d079
1.4/tan(the)

# ╔═╡ 30bd5410-cc83-43fe-a343-0dc7ac12d2c8
begin
	x_range = -5:0.001:5
	y_range = -0.5:0.001:0.5
	
	grid_points = collect(Iterators.product(x_range, y_range))

	for position
	
end

# ╔═╡ 29bea1aa-62c3-4b8f-bf92-6f6e0792e98d
begin 
df = CSV.File("/Users/pavelzhelnin/Documents/physics/TAMBO/notebooks/WhitePaper_EarthMantle.csv")
end 

# ╔═╡ 0052ae3e-6d6d-4ecd-ab21-4c6350313fa5
begin
	result_index = [res for res in df["index"]]
	@bind pass_index Select(result_index)
end

# ╔═╡ c0b7f823-9afa-4bd8-9399-49b7916cec9e
full_particle_viewer(pass_index)

# ╔═╡ 64c69f65-f61e-4ef3-a723-fd097cf08dc2
begin 
inter_azi = atan.(df["ydir"]./df["xdir"])*180/pi

start_index = findall(x -> x >50, inter_azi)

# Shift the elements after start_index by adding 10
inter_azi[start_index] .-= 180

l   = @layout [a{0.75w} b{0.5h}; c{0.6h} ]
plt = histogram(inter_azi,bins=100,layout = l,xlabel="azimuthal angle of passing events",legend=false)
vline!(plt[1],[45],color=:red, linewidth=2, linestyle=:dash)
vline!(plt[1],[-135],color=:red, linewidth=2, linestyle=:dash)

quiver!(plt[2],zeros(size(df["xdir"])),zeros(size(df["xdir"])),quiver =(df["xdir"],df["ydir"]),legend=false,xtickfontsize=3,xlabel="")

print(minimum(inter_azi.-39.5))
quiver!(plt[2],[0],[0],quiver=([0.7071067812],[0.7071067812]),color="red")
quiver!(plt[2],[0],[0],quiver=([-0.7071067812],[-0.7071067812]),color="red")

quiver!(plt[3],df["cc_x"],df["cc_y"],quiver =(df["xdir"],df["ydir"]),legend=false,xlabel="x-y distribution")
end 

# ╔═╡ 9300e1e5-fefd-47c8-af70-e4f3a77cea28
begin 
inter_z = df["inter_z"].*1000

histogram(inter_z, bins=100 ,xlabel="elevation of intercepting events (TAMBO) (m)",labels = "")
vline!([700],color="red",linewidth=2,ylim=(0,100),legend =false,labels = "")
plot!([700,10000],[0,0],fillrange=[100,100],fillalpha = 0.5,color="red",legend = true, labels = "down-going events")
end 

# ╔═╡ 88a23e1e-5980-44a0-a480-cdcfcb6adfdf
begin
	ind = findall(x -> x < 700, df["inter_z"].*1000)
	print(df["index"][ind])
end

# ╔═╡ 9a6a854b-2f35-4b65-a8e8-e49885726797
begin
	f(x,y) = geo(x*units.m,y*units.m)/units.m
	xx = LinRange(400 + geo.box.c1[1]/units.m, geo.box.c2[1]/units.m, 100)
	yy = LinRange(400 + geo.box.c1[2]/units.m, geo.box.c2[2]/units.m, 100)
	
	kwargs = Dict(
	          :st=>:surface,
	          :alpha=>0.9,
	          :c=>cgrad(palette([:skyblue3, :skyblue2, :navajowhite3, :navajowhite3, :goldenrod4, :goldenrod4, :olivedrab, :olivedrab, :green, :green, :green, :green])),
	          :colorbar=>false,
	          :legend=>false,
	         )
	surface(xx, yy, f, #=axis=([], false),=# camera = (azimuth1,90); kwargs...)
	scatter!([df["cc_x"].*1000],[df["cc_y"].*1000],[df["cc_z"].*1000])
end

# ╔═╡ 870f7357-ec9c-480c-a7c5-ac493e3c9574
begin
	surface(xx, yy, f, #=axis=([], false),=# camera = (azimuth1,90); kwargs...)
	scatter!([df["inter_x"].*1000],[df["inter_y"].*1000],[df["inter_z"].*1000])
end

# ╔═╡ 3b27d53b-f856-45c7-b6e0-ac84cb3e4c96
scatter([df["inter_x"].*1000],[df["inter_y"].*1000])

# ╔═╡ 944e3508-d4f5-4e2f-91ad-8630b58186d9
function shower_viewer(index)

	s_i = df["second_index"][index]
	i = df["index"][index]
	
	pos = DataFrame(Dataset("/Users/pavelzhelnin/Documents/physics/TAMBO/notebooks/sim_test_data/sim_test_$(i)_$(s_i)/particles/particles.parquet"))

	position = pos[:, [:x, :y, :z]] 
	
	position.z = (((position.x).*-0.73200153565805759)+((position.y).*0.59897043830965668))/-0.3246662375815868

	position.x.+=(df["inter_x"][index].*1000)
	position.y.+=(df["inter_y"][index].*1000)
	position.z.+=(df["inter_z"][index].*1000)

	CC_pos = SVector(df["cc_x"][index].*1000,df["cc_y"][index].*1000,df["cc_z"][index].*1000)

	taudecay_pos = SVector(df["x"][index].*1000,df["y"][index].*1000,df["z"][index].*1000)

	#indices = findall(x -> -1 <= x <= 1, position.z)
	
	#position = position[indices,:]

	position = Transpose(Matrix(position))

	xx = LinRange(6000 + geo.box.c1[1]/units.m, geo.box.c2[1]/units.m-8000, 100)
	yy = LinRange(6000 + geo.box.c1[2]/units.m, geo.box.c2[2]/units.m-8000, 100)

	kwargs = Dict(
          :st=>:surface,
          :alpha=>0.9,
          :c=>cgrad(palette([:skyblue3, :skyblue2, :navajowhite3, :navajowhite3, :goldenrod4, :goldenrod4, :olivedrab, :olivedrab, :green, :green, :green, :green])),
          :colorbar=>false,
          :legend=>false,
         )
	
	f(x,y) = geo(x*units.m,y*units.m)/units.m

	surface(xx, yy, f, axis=([], false), camera = (-50,20); kwargs...)
	scatter!([position[1,:]],[position[2,:]],[position[3,:]])

	scatter!([df["x"][index].*1000],[df["y"][index].*1000],[df["z"][index].*1000],color="red")

	scatter!([df["cc_x"][index].*1000],[df["cc_y"][index].*1000],[df["cc_z"][index].*1000],color="green")

	track = Tambo.Track(CC_pos,taudecay_pos)

	λλ = LinRange(0, 1, 101)
	λλ = λλ[λλ.<=λcut]
	t = track.(λλ)
	
	plot3d!(
		getindex.(t ,1 ), 
		getindex.(t, 2 ),
		getindex.(t, 3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=false,
		legend=false
	)

	#scatter!([df["inter_x"][index].*1000],[df["inter_y"][index].*1000],[df["inter_z"][index].*1000],color="black")
	
end

# ╔═╡ 8525363e-a06c-4b0a-9547-b5c3f1e88713
print(df["index"])

# ╔═╡ e5f74cda-161c-4851-a9a4-f3498aca6d28
df["index"]

# ╔═╡ 5bf64dc3-6c27-448b-9beb-96a91f1eb0e3
@bind azimuth2 Slider(-180:180; show_value=true, default=15)

# ╔═╡ 1049b43a-56dc-495a-84a9-1fe4ce92d410
@bind zenith2 Slider(0:90; show_value=true, default=15)

# ╔═╡ 5100ba61-0f86-4656-a280-73f0fef59874
begin
	resultr_index = [res for res in df["index"]]
	print(findfirst(x->x==2327,resultr_index))
	index_range = range(1,length(df["index"]))
	@bind pass_pass_index Select(index_range)
end

# ╔═╡ 0f465285-ccc5-43ec-82a6-151049e52b30
shower_viewer(pass_pass_index)

# ╔═╡ a94515ad-bd66-4ee8-9ed3-c936424e793b
function sample_points(n, r, x, y)
    points = []
    
    for i in 1:n
        angle = 2π * rand()  # Randomly sample an angle between 0 and 2π
        distance = sqrt(rand()) * r  # Randomly sample a distance within the radius
        
        # Calculate the coordinates of the sampled point
        point_x = x + distance * cos(angle)
        point_y = y + distance * sin(angle)
        
        push!(points, (point_x, point_y))  # Add the sampled point to the list
    end
    
    return points
end

# ╔═╡ e47dc3ea-9ab5-4f3c-aed4-5c7ffc3c7a00
function test_shower_viewer(azimuth2)

	#s_i = df["second_index"][index]
	#i = df["index"][index]
	
	#pos = DataFrame(Dataset("/Users/pavelzhelnin/Documents/physics/TAMBO/notebooks/sim_test_data/sim_test_$(i)_$(s_i)/particles/particles.parquet"))

	#position = pos[:, [:x, :y, :z]] 
	
	#position.z = (((position.x).*-0.73200153565805759)+((position.y).*0.59897043830965668))/-0.3246662375815868

	#position.x.+=(df["inter_x"][index].*1000)
	#position.y.+=(df["inter_y"][index].*1000)
	#position.z.+=(df["inter_z"][index].*1000)
	#indices = findall(x -> -1 <= x <= 1, position.z)
	
	#position = position[indices,:]

	#position = Transpose(Matrix(position))

	xx = LinRange(9000 + geo.box.c1[1]/units.m, geo.box.c2[1]/units.m-14000, 100)
	yy = LinRange(6000 + geo.box.c1[2]/units.m, geo.box.c2[2]/units.m-24000, 100)

	kwargs = Dict(
          :st=>:surface,
          :alpha=>0.9,
          :c=>cgrad(palette([:skyblue3, :skyblue2, :navajowhite3, :navajowhite3, :goldenrod4, :goldenrod4, :olivedrab, :olivedrab, :green, :green, :green, :green])),
          :colorbar=>true,
          :legend=>false,
         )
	
	f(x,y) = geo(x*units.m,y*units.m)/units.m

	surface(xx, yy, f, axis=([], false), camera = (azimuth2,zenith2); kwargs...)
	
	#scatter!([position[1,:]],[position[2,:]],[position[3,:]])

	CC_pos = SVector(-2000,2000,-100)
	taudecay_pos = SVector(-1000,1000,100)

	scatter!([CC_pos[1]],[CC_pos[2]],[CC_pos[3]],color="red")
	scatter!([taudecay_pos[1]],[taudecay_pos[2]],[taudecay_pos[3]],color="green")
	
	track = Tambo.Track(CC_pos,taudecay_pos)

	λλ = LinRange(0.16, 0.84, 101)
	λλ = λλ[λλ.<=λcut]
	t = track.(λλ)
	
	#print(findfirst(x->Tambo.inside(getindex.(x,1),getindex.(x,2),getindex.(x,3),geo),t))
	
	plot3d!(
		getindex.(t, 1 ), 
		getindex.(t, 2 ),
		getindex.(t, 3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=false,
		legend=false
	)

	λλ = LinRange(-3, -0.16, 101)
	λλ = λλ[λλ.<=λcut]
	bt = track.(λλ)

	plot3d!(
		getindex.(bt, 1 ), 
		getindex.(bt, 2 ),
		getindex.(bt, 3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=false,
		legend=false
	)

	oned_valley(λ) = geo(track(λ).x, track(λ).y)
    root_func(λ) = oned_valley(λ) - track(λ).z
    zeros = Tambo.find_zeros(root_func, 0, 4)
	intersection = track(zeros[1])
	#intersection_track = Tambo.Track(taudecay_pos,intersection)
	#it = intersection_track.(λλ)

	pts = sample_points(100,300,intersection[1],intersection[2])
	x_values = [point[1] for point in pts]
	y_values = [point[2] for point in pts]

	z_values = [geo(x*units.m,y*units.m)/units.m for (x,y) in zip(x_values,y_values)]

	scatter!(x_values,y_values,z_values,marker=:star,color="yellow",markerstrokecolor="yellow",alpha=0.3)

	#Tambo.inside(CC_pos*units.m,geo)

	#=plot3d!(
		getindex.(it ,1 ), 
		getindex.(it, 2 ),
		getindex.(it, 3 ),
		line_z=λλ,
		clim=(0,1),
		lw=2,
		colorbar=false,
		legend=false
	)=#
	
	#scatter!([intersection[1]],[intersection[2]],[intersection[3]],color="black")

	#scatter!([df["inter_x"][index].*1000],[df["inter_y"][index].*1000],[df["inter_z"][index].*1000],color="black")
	
end

# ╔═╡ c10da3bc-f84d-42cd-a655-d797ff444b79
test_shower_viewer(azimuth2)

# ╔═╡ 87f8c4cd-a7df-4b41-89a6-8dafea25b1e9
begin
	anim = @animate for i in -180:179
		test_shower_viewer(i)
	end
	gif(anim, "/Users/pavelzhelnin/Documents/physics/TamboShower.gif", fps = 30)
end

# ╔═╡ 8b62919f-c68d-436a-b23e-905695948b8a
begin 
a = 0.452174
b = -0.366163
c = 0.813304
print(findmax(df["inter_x"]*a + df["inter_y"]*b + df["inter_z"]*c))
print(findmax(df["xdir"]*a + df["ydir"]*b + df["zdir"]*c))
end


# ╔═╡ Cell order:
# ╠═75cc79a1-a15e-4631-a275-3ee3053aaa5d
# ╠═c8b8d3bd-dedf-4333-b0e7-764e7aa2e9cd
# ╠═3ac8cead-b089-4720-91cc-10fef3dbcd5a
# ╠═b08b4b48-4e39-4a3b-858c-47a2666c3b04
# ╠═c0b7f823-9afa-4bd8-9399-49b7916cec9e
# ╠═0052ae3e-6d6d-4ecd-ab21-4c6350313fa5
# ╠═162daaf3-14ba-4099-9308-8f2c7257eede
# ╠═60518038-b4d8-4000-bd9e-d77015217cab
# ╠═bba917d0-f929-4ed7-ad5d-4f08912d7286
# ╠═2889faad-684d-4e93-9fc4-2d8cb07fc7a2
# ╠═0a203549-8d17-43c1-91c5-1e6018544d3f
# ╠═220345e2-83b3-4434-bcb4-3b5f5352b8bf
# ╠═b59160aa-5245-4827-ab0b-c6f73cc05402
# ╠═e64ab6a5-edb6-482d-ab70-0aeeeb425431
# ╠═64c69f65-f61e-4ef3-a723-fd097cf08dc2
# ╠═9300e1e5-fefd-47c8-af70-e4f3a77cea28
# ╠═8c585aae-263a-4ba8-a8bb-4a3afb5f00be
# ╠═88a23e1e-5980-44a0-a480-cdcfcb6adfdf
# ╠═9a6a854b-2f35-4b65-a8e8-e49885726797
# ╠═870f7357-ec9c-480c-a7c5-ac493e3c9574
# ╠═e57cd392-fc79-4359-b461-d795b82ec082
# ╠═25f16065-881f-4d85-9091-70813bb6811b
# ╠═b79157b5-9f07-413c-96cc-cd2905dfd86d
# ╠═fe9c6cea-463b-4474-bcc2-ddeb30b19391
# ╠═58a0ddb1-d1b1-4c87-8e06-5d106df268f4
# ╠═c059d526-77dd-49c7-ab2b-7aa8cf9dbe69
# ╠═9d746940-81ba-4b66-b4f0-e3b843a60cdd
# ╠═1df771f5-50e7-410b-aa85-a5c138c3d7a4
# ╠═3b27d53b-f856-45c7-b6e0-ac84cb3e4c96
# ╠═944e3508-d4f5-4e2f-91ad-8630b58186d9
# ╠═0f465285-ccc5-43ec-82a6-151049e52b30
# ╠═8525363e-a06c-4b0a-9547-b5c3f1e88713
# ╠═ab37cf29-61fe-4f8d-a251-5453e3e4d079
# ╠═e5f74cda-161c-4851-a9a4-f3498aca6d28
# ╠═30bd5410-cc83-43fe-a343-0dc7ac12d2c8
# ╠═29bea1aa-62c3-4b8f-bf92-6f6e0792e98d
# ╠═e47dc3ea-9ab5-4f3c-aed4-5c7ffc3c7a00
# ╠═c10da3bc-f84d-42cd-a655-d797ff444b79
# ╠═5bf64dc3-6c27-448b-9beb-96a91f1eb0e3
# ╠═1049b43a-56dc-495a-84a9-1fe4ce92d410
# ╠═5100ba61-0f86-4656-a280-73f0fef59874
# ╠═a94515ad-bd66-4ee8-9ed3-c936424e793b
# ╠═87f8c4cd-a7df-4b41-89a6-8dafea25b1e9
# ╠═8b62919f-c68d-436a-b23e-905695948b8a
