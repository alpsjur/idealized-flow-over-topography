using Oceananigans
using Statistics
using Dates

function run_model(Nx, Ny, Nz, Δt, stop_time, architecture)
    # grid parameters
    Lx = 2π 
    Ly = 2π 
    Lz = (2π+0.1*2π)*1e-3  

    # misc parameters
    β = 0          # planetary beta
    f = 0          # rotation 
    α = 1e-3       # topographic slope

    # Create grid
    underlying_grid = RectilinearGrid(
        architecture;
        size=(Nx, Ny, Nz), 
        x = (0, Lx),
        y = (0, Ly),
        z = (-Lz, 0), 
        halo = (4, 4, 4),
        topology=(Periodic, Bounded, Bounded)
        )

    # define bathymetry. Slope, depth decreasing towards positive y
    hᵢ(x, y) = -Lz+α*y     

    # create grid with immersed bathymetry 
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(hᵢ))

    # Rotation
    coriolis = BetaPlane(f₀=f, β=β)

    # Turbulence closures 
    κh = 0
    νh = 1e-5   
    κz = 0
    νz = 1e-5
    vertical_closure = ScalarDiffusivity(ν = νz, κ = κz)                 
    horizontal_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)

    # create model
    model = HydrostaticFreeSurfaceModel(; grid,
            free_surface = ImplicitFreeSurface(),
            momentum_advection = WENO(),
            tracer_advection = WENO(),
            closure = (horizontal_closure, vertical_closure),
            coriolis = coriolis
            )

    # set initial random velocity field
    u, v, w = model.velocities

    uᵢ = rand(size(u)[1:2]...)
    vᵢ = rand(size(v)[1:2]...)

    uᵢ .-= mean(uᵢ)
    vᵢ .-= mean(vᵢ)

    uᵢ = repeat(uᵢ, 1, 1, size(u)[3])
    vᵢ = repeat(vᵢ, 1, 1, size(v)[3])

    set!(model, u=uᵢ, v=vᵢ, w=0)

    # create simulations
    simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

    run!(simulation)
    return simulation
end 

function benchmark_model(Nx, Ny, Nz, Δt, stop_time, architecture, initial)
    simulation = run_model(Nx, Ny, Nz, Δt, stop_time, architecture)
    wall_time = simulation.run_wall_time

    # convert variables to string
    snow = string(now())
    shostname = gethostname()
    sarchitecture = string(architecture)[1:3]
    sthreads = string(Threads.nthreads())
    sinitial = (initial)
    sdt = string(Δt)
    ssteps = string(steps)
    sNx = string(Nx)
    sNy = string(Ny)
    sNz = string(Nz)
    stime = string(round(wall_time, digits = 2))
    output = snow*" "*shostname*" "*sarchitecture*" "*sthreads*" "*sinitial*" "*sdt*" "*ssteps*" "*sNx*" "*sNy*" "*sNz*" "*stime


    file =  open("benchmark_results.txt","a")
    write(file, output*"\n")

    return simulation
end    


architecture = GPU()

# simulation parameters
Δt = 0.005
steps = 1000
stop_time = Δt*steps

# run model once to compile model 
run_model(8, 8, 1, Δt, Δt*2, architecture)

# numer of times to run simulatio for each configuration
nsim = 4

for Nh in (8, 16, 32, 64, 128, 256, 512) 
    Nx = Nh
    Ny = Nh
    for Nz in (1, 2, 4, 8, 16, 32)
        initial = true
        for i in 1:nsim
            benchmark_model(Nx, Ny, Nz, Δt, stop_time, architecture, initial)
            initial = false
        end
    end
end
