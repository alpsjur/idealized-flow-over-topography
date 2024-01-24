using Oceananigans
using Oceananigans.Units
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

# bathymetry parameters (see Nummelin & Isachsen, 2024)
W  =   75kilometers
YC =  150kilometers
DS = 2000meters
DB =  250meters
σ  =   10meters


# simulation parameters 
Δt = 120seconds
stop_time = 365days
save_fields_interval = 1days

# forcing parameters
τ = 0.05/1000    # N m-2 / kg m-2   kinematic forcing

# bottom friction?
r = 3e-4            # Bottom linear drag coefficient m/s

# stratification parameters TODO sjekk dette
Tmax = 15
decay = 50meter 

# Create grid
underlying_grid = RectilinearGrid(
                                CPU();
                                size=(Nx, Ny, Nz), 
                                x = (0, Lx),
                                y = (0, Ly),
                                z = (-Lz,0), 
                                halo = (4, 4, 4),
                                topology=(Periodic, Bounded, Bounded)
                                )

# define bathymetry
function hᵢ(x, y)
    if y < (YC + W)                # shelf
        h = -DB - 0.5*DS*(1+tanh.(π*(y-YC)/W))
    elseif Ly - y < (YC + W)       # slope
        h = -DB - 0.5*DS*(1+tanh.(π*(Ly-y-YC)/W))
    else                           # central basin
        h = -DB - DS
    end

    # add random noise
    #h .+= randn((length(x),length(y)))*σ
  
    return h
end


# create grid with immersed bathymetry 
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(hᵢ))

println(grid)

# Turbulence closures 
κh = 0
νh = 1e-5   
κz = 0
νz = 1e-5
vertical_closure = ScalarDiffusivity(ν = νz, κ = κz)                 
horizontal_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)

# Rotation
coriolis = FPlane(1e-4)


# spesify bottom drag 
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
                                    coriolis = coriolis
                                    )

    
# define initial temperature TODO update
initial_temperature(x, y, z) = Tmax*exp(-decay*z)

# set initial temperature
set!(model, T=initial_temperature)                    
                   
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
          
                                                                
# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 5.2f, wall time: % 15s, max |u|: % 5.3f, max |v|: % 5.3f, max |w|: % 5.3f\n",
                        sim.model.clock.iteration,
                        #prettytime(sim.model.clock.time),
                        sim.model.clock.time,
                        prettytime(1e-9 * (time_ns() - start_time)),
                        maximum(abs, sim.model.velocities.u),
                        maximum(abs, sim.model.velocities.v),
                        maximum(abs, sim.model.velocities.w),
                        )

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

# write output to file
filename = "BLOM_channel"
datapath = "channel/data/"

u, v, w = model.velocities

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w),
                                                      schedule = TimeInterval(save_fields_interval),
                                                      filename = datapath*filename*".jld2",
                                                      overwrite_existing = true
                                                      )


# action!
run!(simulation)