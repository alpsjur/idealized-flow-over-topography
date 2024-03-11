include("channel_setup.jl")

# TODO 
# save bottom drag 

# Overwrite variables from channel_setup.jl
Nx = 20
Lx = dx*Nx
stop_time = 100days

# Create grid
underlying_grid = RectilinearGrid(
        architecture;
        size=(Nx, Ny, Nz), 
        x = (0, Lx),
        y = (0, Ly),
        #y = y_faces,
        z = z_faces,
        halo = (4, 4, 4),
        topology=(Periodic, Bounded, Bounded)
)
                              
# create grid with immersed bathymetry 
grid = ImmersedBoundaryGrid(underlying_grid, 
                            #GridFittedBottom(hᵢ)
                            PartialCellBottom(hᵢ, minimum_fractional_cell_height=0.1)
                            )


"""
# visualize vertical grid spacing
fig = Figure()
ax = Axis(fig[1, 1], ylabel = "Depth (m)", xlabel = "Vertical spacing (m)")

lines!(ax, zspacings(underlying_grid, Center()), znodes(underlying_grid, Center()))
scatter!(ax, zspacings(underlying_grid, Center()), znodes(underlying_grid, Center()))

ax = Axis(fig[1, 2], xlabel = "k index")
ks = collect(1:51)
zs = z_faces.(ks)

scatter!(ax, ks, zs)
save(figurepath*"vertical_grid_spacing.png", fig)

# visualize y grid spacing
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "y (m)", ylabel = "y spacing (m)")

lines!(ax, ynodes(underlying_grid, Center()), yspacings(underlying_grid, Center()))
scatter!(ax, ynodes(underlying_grid, Center()), yspacings(underlying_grid, Center()))
save(figurepath*"y_grid_spacing.png", fig)

# visualize bathymetry
x, y, z = nodes(grid, (Center(), Center(), Center()))

bath = grid.immersed_boundary.bottom_height.data[1:Nx,1:Ny]

fig, ax, hm = heatmap(x*1e-3, y*1e-3, bath,
                      colormap=:deepsea,
                      axis = (xlabel = "x [km]",
                              ylabel = "y [km]",
                              title = "Bathymetry",
                              titlesize = 24))

Colorbar(fig[1, 2], hm, label = "depth [m]")

current_figure() # hide
save(figurepath*"channel_bathymetry.png", fig)  
"""

# create model
model = HydrostaticFreeSurfaceModel(; 
        grid,
        #grid=underlying_grid,
        boundary_conditions=(u=u_bc, v=v_bc),
        free_surface = ImplicitFreeSurface(),
        momentum_advection = WENO(),
        tracer_advection = WENO(),
        closure = closure,
        coriolis = coriolis,
        buoyancy = BuoyancyTracer(),
        tracers = :b,
)

println(model)

# set initial density profile
set!(model, b=initial_buoyancy)             

"""
# plot initial profile
fig = Figure()
axis = Axis(fig[1,1], xlabel = "Buoyancy", ylabel = "z")

z = znodes(model.tracers.b)
b = interior(model.tracers.b, 104, 256, :)

lines!(axis, b, z)
save(figurepath*"initial_buoyancy.png", fig)
"""                      
                                                                
# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# test time step wizard
#wizard = TimeStepWizard(cfl=0.2)
#simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))


#define diagnostics  (Which to save?)
#define diagnostics 
include("diagnostics.jl")
#output_attributes=output_attributes,


# logging simulation progress
simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))

# write output to file
filename = "andre_closure"
datapath = "channel/data/run5/"

Uh = Average(u, dims=1)
Vh = Average(v, dims=1)
Wh = Average(w, dims=1)
Bh = Average(b, dims=1)

Uv = Average(u, dims=3)
Vv = Average(v, dims=3)
Wv = Average(w, dims=3)
Bv = Average(b, dims=3)

H = Average(η′, dims=(1,3))

"""
simulation.output_writers[:fields] = JLD2OutputWriter(
        model, (; 
                u, v, w,
                #uu, vv, uv,
                η, 
                p, 
                #p_b,
                b,
                #ub, vb, wb,  
        ),
        schedule = AveragedTimeInterval(
                save_fields_interval, 
                window=average_window
        ),
        filename = datapath*filename*".jld2",
        overwrite_existing = true,
        with_halos = true,                           # for computation of derivatives at boundaries
        init = init_save_some_metadata!
)
"""

simulation.output_writers[:averages] = JLD2OutputWriter(
        model, (; 
                Uh, Vh, Wh,
                Bh, 
                H, 
                Uv, Vv, Wv,
                Bv,
        ),
        schedule = AveragedTimeInterval(
                save_fields_interval, 
                window=average_window
        ),
        filename = datapath*filename*"_average"*".jld2",
        overwrite_existing = true,
        with_halos = true,                           # for computation of derivatives at boundaries
        init = init_save_some_metadata!
)

"""
# Remove netCDF file if it already exists 
if isfile(datapath*filename*".nc")
        rm(datapath*filename*".nc")
end
    

# Writer for netCDF file format
simulation.output_writers[:netCDF] = NetCDFOutputWriter(
        model, 
        Dict(
                "u" => u,
                "v" => v,
                "w" => w,
                "b" => b,
                "p" => p,
                "eta" => η,
        ),
        output_attributes=output_attributes,
        schedule = AveragedTimeInterval(
                save_fields_interval, 
                window=average_window
        ),
        filename = datapath*filename*".nc",
        overwrite_existing = true,
        with_halos = true,                     # for computation of derivatives at boundaries. 
)

# Writer for netCDF file format
simulation.output_writers[:grid] = NetCDFOutputWriter(
        model, 
        Dict(
                "h" => h,
        ),
        dimensions=dims,
        output_attributes=output_attributes,
        schedule = SpecifiedTimes([0]),
        filename = datapath*filename*"_bottom_height.nc",
        overwrite_existing = true,
        with_halos = true,                     
)

"""

# action!
run!(simulation)
