include("channel_setup.jl")



# Create grid
underlying_grid = RectilinearGrid(
  architecture;
  size=(Ny, Nz), 
  y = (0, Ly),
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
  grid,
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



"""
### Code from https://github.com/CliMA/Oceananigans.jl/discussions/3423
fcf = (Face(), Center(), Face())


get_qᵂ(i, j, k, ibg, args...)
get_qᴱ(i, j, k, ibg, args...)
get_qˢ(i, j, k, ibg, args...)
get_qᴺ(i, j, k, ibg, args...)
get_qᴮ(i, j, k, ibg, args...)
get_qᵀ(i, j, k, ibg, args...)

 @inline function immersed_flux_divergence(i, j, k, ibg::GFIBG, bc, loc, c, closure, K, id, clock, fields) 
     return div(i, j, k, ibg, loc, get_qᵂ(...), get_qᴱ(...), get_qˢ(...), get_qᴺ(...), get_qᴮ(...), get_qᵀ(...)) 
 end 

using Oceananigans.ImmersedBoundaries: conditional_flux, bottom_ib_flux
using Oceananigans.Operators: index_left
using Oceananigans.BoundaryConditions: flip
@inline function conditional_bottom_ib_flux(i, j, k, ibg, bc, loc, c, closure, K, id, clock, fields)
    q̃ᴮ = bottom_ib_flux(i, j, k, ibg, bc.bottom, loc, c, closure, K, id, clock, fields)

    iᵂ, jˢ, kᴮ = map(index_left,  (i, j, k), loc) # Broadcast instead of map causes inference failure
    LX, LY, LZ = loc
    return conditional_flux(i, j, kᴮ, ibg, LX, LY, flip(LZ), q̃ᴮ, zero(eltype(ibg)))
end

τᵤᶻ_ib = Field(KernelFunctionOperation{Face, Center, Face}(conditional_bottom_ib_flux, grid,
                                                          u.boundary_conditions.immersed, fcf, u, model.closure,
                                                          model.diffusivity_fields, nothing, model.clock, fields(model)))
compute!(τᵤᶻ_ib)

using Oceananigans.Grids: boundary_node
boundary_node_ccf(i, j, k, grid) = boundary_node(i, j, k, grid, Center(), Center(), Face())
boundary_node_ccf_op = KernelFunctionOperation{Center, Center, Face}(boundary_node_ccf, grid)
bn = Field(boundary_node_ccf_op)

τb = Average(τᵤᶻ_ib, dims=3, condition=boundary_node_ccf_op)


###
"""

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf(
  "i: % 6d, sim time: % 15s, wall time: % 15s, max |u|: % 5.3f, max |v|: % 5.3f, max |w|: % 5.3f, max |η|: % 5.3f, next Δt: %s\n",
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
filename = "2D_BLOM_channel"
datapath = "channel/data/"

simulation.output_writers[:fields] = JLD2OutputWriter(
  model, (; 
    u, v, w,
    uu, vv, uv,
    η, p, b,
    ub, vb, wb,  
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