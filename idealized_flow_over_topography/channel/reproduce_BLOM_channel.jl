using Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, PartialCellBottom
using Random
using Printf

# grid parameters 
Lx =  416kilometers
Ly = 1024kilometers
Lz = 2250meters

dx = 2kilometers
dy = 2kilometers


Nx = Int(Lx/dx)
Ny = Int(Ly/dy)
Nz = 51


# TODO fikse z-spacing med økt oppløsning ved overflaten (og ved bunnen?)
refinement = 1.2 # controls spacing near surface (higher means finer spaced)
stretching = 12  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)




# bathymetry parameters (see Nummelin & Isachsen, 2024)
W  =   75kilometers
YC =  150kilometers
DS = 2000meters
DB =  250meters
σ  =   10meters              # std of random noise on topography


# simulation parameters 
Δt = 120second#120seconds
stop_time = Δt*10#2days
save_fields_interval = 1hour

# forcing parameters
τ = -0.05/1000               # N m-2 / kg m-2  kinematic forcing

# bottom friction
r = 3e-4                    # Bottom linear drag coefficient m/s

# stratification parameters 
# svak restoring til denne profilen pga eksplisitt horisontalmiksing + numerisk diffusjon
Tmax = 15
bmax = 
decay = 100meter 

# Create grid
underlying_grid = RectilinearGrid(
                                CPU();
                                size=(Nx, Ny, Nz), 
                                x = (0, Lx),
                                y = (0, Ly),
                                z = (-Lz,0), 
                                #z = z_faces,
                                halo = (4, 4, 4),
                                topology=(Periodic, Bounded, Bounded)
                                )
                              

# define bathymetry
function hᵢ(x, y)
    if y < (YC + W)                # shelf
        h =  -DB - 0.5*DS*(1+tanh.(π*(y-YC)/W))
    elseif Ly - y < (YC + W)       # slope
        h = -DB - 0.5*DS*(1+tanh.(π*(Ly-y-YC)/W))
    else                           # central basin
        h =  -DB - DS
    end
    # add random noise
    h += randn()*σ
    return h
end


# create grid with immersed bathymetry 
grid = ImmersedBoundaryGrid(underlying_grid, 
                            #GridFittedBottom(hᵢ)
                            PartialCellBottom(hᵢ)
                            )
#grid = underlying_grid
println(grid)

# visualize bathymetry
using CairoMakie

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
save("channel_bathymetry.png", fig)  


# Turbulence closures 
κh = 0
νh = 0   
κz = 1e-5
νz = 1e-5
vertical_closure = ScalarDiffusivity(ν = νz, κ = κz)                 
horizontal_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)

# Rotation
coriolis = FPlane(1e-4)


# spesify bottom drag TODO Vurder kvadratisk, r->R
drag_u(x, y, t, u, v, r) = -r*u 
drag_v(x, y, t, u, v, r) = -r*v 

drag_u_bc = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=r)
drag_v_bc = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=r)

# spesify drag on immersed boundary, note different arguments from bottom drag above
immersed_drag_u(x, y, z, t, u, v, r) = -r*u 
immersed_drag_v(x, y, z, t, u, v, r) = -r*v 

immersed_drag_u_bc = FluxBoundaryCondition(immersed_drag_u, field_dependencies=(:u, :v), parameters=r)
immersed_drag_v_bc = FluxBoundaryCondition(immersed_drag_v, field_dependencies=(:u, :v), parameters=r)

immersed_u_bc = ImmersedBoundaryCondition(bottom = immersed_drag_u_bc)
immersed_v_bc = ImmersedBoundaryCondition(bottom = immersed_drag_v_bc)

# spesify surface forcing
τx(x, y, t, τ) = τ
τy(x, y, t, τ) = 0

τx_bc = FluxBoundaryCondition(τx, parameters=τ)
τy_bc = FluxBoundaryCondition(τy, parameters=τ)

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
                                    #buoyancy = BuoyancyTracer(),
                                    #tracers = :b,
                                    )

println(model)

# TODO hvilken equation of state er i Aleksis simuleringer? 
    
# define initial temperature TODO update
initial_temperature(x, y, z) = Tmax.*exp.(z./decay)
initial_buoyancy(x, y, z) = bmax.*exp.(z./decay)

# set initial density profile
set!(model, T=initial_temperature, S=34)  
#set!(model, b=initial_buoyancy)              

# plot initial profile
fig = Figure()
axis = Axis(fig[1,1], xlabel = "Temperature (ᵒC)", ylabel = "z")

z = znodes(model.tracers.T)
T = interior(model.tracers.T, 1, 1, :)

lines!(axis, T, z)
save("initial_temperature.png", fig)
                      
                                                                
# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# test time step wizard
wizard = TimeStepWizard(cfl=0.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 5.2f, wall time: % 15s, max |u|: % 5.3f, max |v|: % 5.3f, max |w|: % 5.3f, next Δt: %s\n",
                        sim.model.clock.iteration,
                        #prettytime(sim.model.clock.time),
                        sim.model.clock.time,
                        prettytime(1e-9 * (time_ns() - start_time)),
                        maximum(abs, sim.model.velocities.u),
                        maximum(abs, sim.model.velocities.v),
                        maximum(abs, sim.model.velocities.w),
                        prettytime(sim.Δt),
                        )

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

# write output to file
filename = "BLOM_channel"
datapath = "channel/data/"

# TODO kan tidsmidlede variabler lagres?
u, v, w = model.velocities

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w),
                                                      schedule = TimeInterval(save_fields_interval),
                                                      filename = datapath*filename*".jld2",
                                                      overwrite_existing = true
                                                      )


# action!
run!(simulation)