include("channel_setup.jl")
figurepath = "channel/figures/"

Nx = 10
Ny = 40

Lx = dx*Nx 
Ly = dy*Ny

stop_time = 100days 

function hᵢ(x, y)
    if y < Ly/2 
        return -Lz
    else 
        return -Lz + (y-Ly/2)*0.01
    end
end

# Create grid
underlying_grid = RectilinearGrid(
        architecture;
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



"""
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
save(figurepath*"test_bottom_drag_bathymetry.png", fig)  
"""

# create model
model = HydrostaticFreeSurfaceModel(; 
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

println(model)

# set initial density profile
#set!(model, b=initial_buoyancy)  
  

# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# test time step wizard
#wizard = TimeStepWizard(cfl=0.2)
#simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))


#define diagnostics 
include("diagnostics.jl")

# logging simulation progress
simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))


# write output to file
filename = "test_bottom_drag_output"
datapath = "channel/data/"


simulation.output_writers[:fields] = JLD2OutputWriter(
        model, (; 
        u, v, w,
        u_bc_op, v_bc_op, u_im_bc_op, v_im_bc_op,
        u_b, v_b,
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