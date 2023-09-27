using Oceananigans
using Statistics
using CairoMakie
using Printf

# grid parameters
aspect_ratio = 0.01

Nx = 128
Ny = 128
Nz = 1

Lx = 2π
Ly = 2π
Lz = Lx*aspect_ratio

# misc parameters
β = 0          # planetary beta
ν = 1e-5       # diffusivity

# simulation parameters
Δt = 0.05 
stop_time = 600

# create grid
grid = RectilinearGrid(size=(Nx, Ny, Nz), 
                       extent=(Lx, Ly, Lz), 
                       topology=(Periodic, Bounded, Bounded)
                       )

# model parameters 
coriolis = BetaPlane(f₀=1, β=β)
closure = ScalarDiffusivity(ν=ν)

# create model
model = HydrostaticFreeSurfaceModel(; grid,
                                    closure = closure,
                                    #coriolis = coriolis
                                    )


# set initial random velocity field
u, v, w = model.velocities

uᵢ = rand(size(u)...)
vᵢ = rand(size(v)...)

uᵢ .-= mean(uᵢ)
vᵢ .-= mean(vᵢ)

set!(model, u=uᵢ, v=vᵢ, w=0)

# create simulations
simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 12s, wall time: % 15s, max |u|: % 5s\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(1e-9 * (time_ns() - start_time)),
                        maximum(abs, sim.model.velocities.u),
                        )

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# write output to file
filename = "2Dturbulence_3Ddomain"
datapath = "idealized_flow_over_topography/beta_plane_turbulence/data/"

u, v, w = model.velocities

ω = ∂x(v) - ∂y(u)
Ω = Average(ω, dims=3)

s = sqrt(u^2+v^2)
S = Average(s, dims=3)

U = Average(u, dims=(1, 3))
PV = Average(ω, dims=(1, 3))



simulation.output_writers[:fields] = JLD2OutputWriter(model, (; Ω, S, U, PV),
                                                      schedule = TimeInterval(0.3),
                                                      filename = datapath*filename*".jld2",
                                                      overwrite_existing = true
                                                      )


# action!
run!(simulation)

