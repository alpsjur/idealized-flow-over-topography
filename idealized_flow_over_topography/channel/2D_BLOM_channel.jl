include("channel_setup.jl")

# Create grid
underlying_grid = RectilinearGrid(
                                CPU();
                                size=(Ny, Nz), 
                                y = (0, Ly),
                                z = z_faces,
                                halo = (4, 4),
                                topology=(Flat, Bounded, Bounded)
                                )
                              

# create grid with immersed bathymetry 
# PartialCell dosen't seem to work with 2D model
grid = ImmersedBoundaryGrid(underlying_grid, 
                            GridFittedBottom(hᵢ)
                            #PartialCellBottom(hᵢ)
                            )

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


# create model. 2D model must be nonhydrostatic?
model =  HydrostaticFreeSurfaceModel(; grid,
                                    boundary_conditions=(u=u_bc, v=v_bc),
                                    free_surface = ImplicitFreeSurface(),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    closure = (horizontal_closure, vertical_closure),
                                    coriolis = coriolis,

                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    )


set!(model, b=initial_buoyancy)  

# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# test time step wizard
wizard = TimeStepWizard(cfl=0.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
                   
#define diagnostics  (Which to save?)
u, v, w = model.velocities
p = model.pressure.pHY′      # see here: https://github.com/CliMA/Oceananigans.jl/discussions/3157
η = model.free_surface.η
b = model.tracers.b

uu = u*u 
vv = v*v
uv = u*v
ub = u*b
vb = v*b 
wb = w*b

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 15s, wall time: % 15s, max |u|: % 5.3f, max |v|: % 5.3f, max |w|: % 5.3f, max |η|: % 5.3f, next Δt: %s\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        #sim.model.clock.time,
                        prettytime(1e-9 * (time_ns() - start_time)),
                        maximum(abs, u),
                        maximum(abs, v),
                        maximum(abs, w),
                        maximum(abs, η),
                        prettytime(sim.Δt),
                        )

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# write output to file
filename = "2D_BLOM_channel_test"
datapath = "channel/data/"


simulation.output_writers[:fields] = JLD2OutputWriter(model, (; 
                                                        u, v, w,
                                                        uu, vv, uv,
                                                        η, p, b,
                                                        ub, vb, wb,  
                                                        ),
                                                      schedule = AveragedTimeInterval(save_fields_interval, 
                                                                                      window=average_window),
                                                      filename = datapath*filename*".jld2",
                                                      overwrite_existing = true,
                                                      with_halos = true                           # for computation of derivatives at boundaries
                                                      )


# action!
run!(simulation)