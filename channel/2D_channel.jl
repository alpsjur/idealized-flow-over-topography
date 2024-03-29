include("channel_setup.jl")

stop_time = 50days

# Create grid
underlying_grid = RectilinearGrid(
  architecture;
  size=(Ny, Nz), 
  y = (0, Ly),
  #y = y_faces,
  z = z_faces,
  halo = (4, 4),
  topology=(Flat, Bounded, Bounded)
)
                              

# create grid with immersed bathymetry 
# PartialCell dosen't seem to work with 2D model
grid = ImmersedBoundaryGrid(
  underlying_grid, 
  GridFittedBottom(hᵢ)
  #PartialCellBottom(hᵢ)
)



# create model
model =  HydrostaticFreeSurfaceModel(; 
  #grid,
  grid=underlying_grid,
  boundary_conditions=(u=u_bc, v=v_bc),
  free_surface = ImplicitFreeSurface(),
  momentum_advection = WENO(),
  tracer_advection = WENO(),
  closure = closure,
  coriolis = coriolis,
  buoyancy = BuoyancyTracer(),
  tracers = :b,
)

#set!(model, b=initial_buoyancy)  

println(model)

"""
# plot initial profile
figurepath = "channel/figures/"
fig = Figure()
axis = Axis(fig[1,1], xlabel = "Buoyancy [m/s²]", ylabel = "z [m]")


z = znodes(model.tracers.b)
b = interior(model.tracers.b, 1, 256, :)


lines!(axis, b, z)
save(figurepath*"initial_buoyancy.png", fig)
"""

# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# test time step wizard
#wizard = TimeStepWizard(cfl=0.2)
#simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
                   
#define diagnostics 
include("diagnostics.jl")


# logging simulation progress
simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))

# write output to file
filename = "2D_nostrat_flat_scalar"
datapath = "channel/data/run2/"

simulation.output_writers[:fields] = JLD2OutputWriter(
  model, (; 
    u, v, w,
    uu, vv, uv,
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


# action!
run!(simulation)