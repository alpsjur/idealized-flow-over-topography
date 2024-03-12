using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.Models.HydrostaticFreeSurfaceModels: FFTImplicitFreeSurfaceSolver
using Printf

stop_time = 10days
save_interval = 12hours

#Nx, Nz = 100, 100
#Lx, Lz = 100kilometers, 100meters
Nx, Nz = 150, 150
Lx, Lz = 150kilometers, 150meters
h0 = 30#10
slope = 1e-3

Ny = 10
Ly = Ny*Lx/Nx

Δt = 3minutes
ν = 1e-2
f0 = -5.46e-05
tauy = -1e-4
wind_time_ramp = 0.25days

underlying_grid = RectilinearGrid(GPU(),
                                  topology = (Bounded, Periodic, Bounded),
                                  size = (Nx, Ny, Nz),
                                  x = (0, Lx),
                                  y = (0, Ly),
                                  z = (-Lz, 0))

# An inclined plane beach shelf model.
bathymetry(x, y, z) = z < - (h0 + slope*x)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(bathymetry))
#grid = underlying_grid

wind_stress(x, y, t) = - tauy * tanh(t/wind_time_ramp)
#wind_stress(x, y, t) = 0
no_slip = ValueBoundaryCondition(0)
no_flux = FluxBoundaryCondition(0)

# ####
# #### Linear bottom drag.
# ####
# filename = "slope_2D-linear.nc"
#
# Δz_bottom = grid.Δzᵃᵃᶜ[1]
# r = 1 / 2days * Δz_bottom
# @inline u_bottom_drag(x, y, t, u, v, p) = @inbounds - p.r * u
# @inline v_bottom_drag(x, y, t, u, v, p) = @inbounds - p.r * v
#
# u_immersed_drag_bc = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(u_immersed_drag, field_dependencies=(:u, :v), parameters = (; r)))
# v_immersed_drag_bc = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(v_immersed_drag, field_dependencies=(:u, :v), parameters = (; r)))
# u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, field_dependencies=(:u, :v), parameters = (; r))
# v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, field_dependencies=(:u, :v), parameters = (; r))

####
#### Quadratic bottom friction (based on https://clima.github.io/OceananigansDocumentation/stable/generated/tilted_bottom_boundary_layer/).
####
filename = "slope_3D-quadratic.nc"
#z₀ = 0.1 # m (roughness length)
#κ = 0.4 # von Karman constant
#z₁ = underlying_grid.Δzᵃᵃᶜ[1] # Closest grid center to the bottom
#cᴰ = (κ / log(z₁ / z₀))^2 # Drag coefficient
cᴰ = 3e-3 # Canonical shelf value?

@inline u_immersed_drag(x, y,  z, t, u, v, cᴰ) = @inbounds - cᴰ * u * (u^2 + v^2)^0.5
@inline v_immersed_drag(x, y, z, t, u, v, cᴰ) = @inbounds - cᴰ * v * (u^2 + v^2)^0.5
@inline u_bottom_drag(x, y, t, u, v, cᴰ) = @inbounds - cᴰ * u * (u^2 + v^2)^0.5
@inline v_bottom_drag(x, y, t, u, v, cᴰ) = @inbounds - cᴰ * v * (u^2 + v^2)^0.5

u_immersed_drag_bc = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(u_immersed_drag, field_dependencies=(:u, :v), parameters =cᴰ))
v_immersed_drag_bc = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(v_immersed_drag, field_dependencies=(:u, :v), parameters =cᴰ))
u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, field_dependencies=(:u, :v), parameters = cᴰ)
v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, field_dependencies=(:u, :v), parameters = cᴰ)

#---

bcs = (;
   u = FieldBoundaryConditions(bottom=u_bottom_drag_bc, immersed=u_immersed_drag_bc),
   v = FieldBoundaryConditions(top=FluxBoundaryCondition(wind_stress), bottom=v_bottom_drag_bc, immersed=v_immersed_drag_bc)
)

model = HydrostaticFreeSurfaceModel(grid = grid, 
                                    boundary_conditions=bcs,
                                    momentum_advection = CenteredSecondOrder(),
                                    free_surface = ImplicitFreeSurface(gravitational_acceleration=9.81),
                                    closure = ScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=ν, κ=1e-6),
                                    tracers = :b,
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(f=f0),
                                    )

# Linear stratification
dbdz = 3.2e-4 # = 3.3 kg/m3 top-bottom densitydifference over 100 m.
set!(model, b = (x, y, z) -> dbdz * z)
#set!(model, b = (x, y, z) -> dbdz * z + 1e-6 * x)
#set!(model, b = (x, y, z) -> 4 * z)
#set!(model, b = (x, y, z) -> 4 * z - 1e-8 * x)

progress_message(s) = @info @sprintf("[%.2f%%], time: %.2f days, max|u|: %.3f, max|v|: %.3f",
                            100 * s.model.clock.time / s.stop_time,
                            s.model.clock.time/86400, maximum(abs, model.velocities.u), maximum(abs, model.velocities.v))

gravity_wave_speed = sqrt(model.free_surface.gravitational_acceleration * grid.Lz)
#Δt = 0.1 * minimum(grid.Δxᶜᵃᵃ) / gravity_wave_speed
#Δt = 0.25 * minimum(grid.Δxᶜᵃᵃ) / gravity_wave_speed
simulation = Simulation(model, Δt = Δt, stop_time = stop_time)
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))
# Build the output_writer for the two-dimensional fields to be output.
output = Dict("v" => model.velocities.v, "b" => model.tracers.b)
simulation.output_writers[:fields] = NetCDFOutputWriter(model, output,
                                                        filename = filename,
                                                        schedule = TimeInterval(save_interval),
                                                        overwrite_existing = true)

run!(simulation)
