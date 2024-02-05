include("channel_setup.jl")

# TODO 
# z-spacing
# set averaging time 
# set friction parameter
# save bottom drag 
# save bottom pressure


# Create grid
underlying_grid = RectilinearGrid(
                                CPU();
                                size=(Nx, Ny, Nz), 
                                x = (0, Lx),
                                y = (0, Ly),
                                z = z_faces,
                                halo = (4, 4, 4),
                                topology=(Periodic, Bounded, Bounded)
                                )
                              


# create grid with immersed bathymetry 
grid = ImmersedBoundaryGrid(underlying_grid, 
                            #GridFittedBottom(hᵢ)
                            PartialCellBottom(hᵢ)
                            )


# visualize vertical grid spacing
fig = Figure()
ax = Axis(fig[1, 1], ylabel = "Depth (m)", xlabel = "Vertical spacing (m)")

lines!(ax, zspacings(underlying_grid, Center()), znodes(underlying_grid, Center()))
scatter!(ax, zspacings(underlying_grid, Center()), znodes(underlying_grid, Center()))

ax = Axis(fig[1, 2], xlabel = "k index")
ks = collect(1:51)
zs = z_faces.(k)

scatter!(ax, ks, zs)
save("figures/vertical_grid_spacing.png", fig)


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
save("figures/channel_bathymetry.png", fig)  



drag_u_bc = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=R)
drag_v_bc = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=R)

immersed_drag_u_bc = FluxBoundaryCondition(immersed_drag_u, field_dependencies=(:u, :v), parameters=R)
immersed_drag_v_bc = FluxBoundaryCondition(immersed_drag_v, field_dependencies=(:u, :v), parameters=R)

immersed_u_bc = ImmersedBoundaryCondition(bottom = immersed_drag_u_bc)
immersed_v_bc = ImmersedBoundaryCondition(bottom = immersed_drag_v_bc)


# create boundary conditions
u_bc = FieldBoundaryConditions(
                               bottom=drag_u_bc, 
                               immersed=immersed_u_bc, 
                               top=τx_bc
                               )
v_bc = FieldBoundaryConditions(
                               bottom=drag_v_bc, 
                               immersed=immersed_v_bc, 
                               top=τy_bc
                               )


# create model
model = HydrostaticFreeSurfaceModel(; grid,
                                    boundary_conditions=(u=u_bc, v=v_bc),
                                    free_surface = ImplicitFreeSurface(),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    closure = (horizontal_closure, vertical_closure),
                                    coriolis = coriolis,
                                    
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    )

println(model)

# set initial density profile
#set!(model, T=initial_temperature, S=34)
set!(model, b=initial_buoyancy)              

# plot initial profile
fig = Figure()
axis = Axis(fig[1,1], xlabel = "Buoyancy", ylabel = "z")



z = znodes(model.tracers.b)
b = interior(model.tracers.b, 104, 256, :)


lines!(axis, b, z)
save("figures/initial_buoyancy.png", fig)
                      
                                                                
# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# test time step wizard
wizard = TimeStepWizard(cfl=0.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))


#define diagnostics  (Which to save?)
u, v, w = model.velocities
p = model.pressure      # see here: https://github.com/CliMA/Oceananigans.jl/discussions/3157
η = model.free_surface.η
b = model.tracers

uu = u*u 
vv = v*v
uv = u*v
ub = u*b
vb = v*b 
wb = w*b

# τbx
# τby


# logging simulation progress
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 5.2f, wall time: % 15s, max |u|: % 5.3f, max |w|: % 5.3f, max |η|: % 5.3f, next Δt: %s\n",
                        sim.model.clock.iteration,
                        #prettytime(sim.model.clock.time),
                        sim.model.clock.time,
                        prettytime(1e-9 * (time_ns() - start_time)),
                        maximum(abs, u),
                        maximum(abs, w),
                        maximum(abs, η),
                        prettytime(sim.Δt),
                        )

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

# write output to file
filename = "BLOM_channel"
datapath = "channel/data/"


simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w),
                                                      schedule = AveragedTimeInterval(save_fields_interval, 
                                                                                      window=average_window),
                                                      filename = datapath*filename*".jld2",
                                                      overwrite_existing = true,
                                                      with_halos = true                           # for computation of derivatives at boundaries
                                                      )


# action!
run!(simulation)