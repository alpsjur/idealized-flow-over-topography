using Oceananigans
using Oceananigans.Units
using Statistics
using CairoMakie
using Printf

Nx = 128
Ny = 128
Nz = 8

Lx = 2π 
Ly = 2π 
Lz = 2π*1e-3    # depth [m]

# misc parameters
β = 0          # planetary beta
f = 0         # rotation 

# simulation parameters
Δt = 0.005
stop_time = 500
save_fields_interval = 0.3

# Create grid
grid = RectilinearGrid(size=(Nx, Ny, Nz), 
                       extent=(Lx, Ly, Lz), 
                       topology=(Periodic, Bounded, Bounded)
                       )

# Rotation
coriolis = BetaPlane(f₀=f, β=β)
#coriolis = FPlane(f)

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

# initial field from 2D simulation
#ti = 40

#uᵢ = FieldTimeSeries(datapath*"2D_turbulence_β=0.jld2", "u").data[:,:,1,ti]
#vᵢ = FieldTimeSeries(datapath*"2D_turbulence_β=0.jld2", "v").data[:,:,1,ti]

#uᵢ = repeat(uᵢ, 1, 1, size(u)[3])
#vᵢ = repeat(vᵢ, 1, 1, size(v)[3])

set!(model, u=uᵢ, v=vᵢ, w=0)

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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# write output to file
filename = "2Dturbulence_3Ddomain_β=$β"
datapath = "idealized_flow_over_topography/beta_plane_turbulence/data/"
animationpath = "idealized_flow_over_topography/beta_plane_turbulence/animations/"

u, v, w = model.velocities

ω = ∂x(v) - ∂y(u)
Ω = Average(ω, dims=3)

s = sqrt(u^2+v^2)
S = Average(s, dims=3)

U = Average(u, dims=(1, 3))
PV = Average(ω, dims=(1, 3))



simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, Ω, S, U, PV),
                                                      schedule = TimeInterval(save_fields_interval),
                                                      filename = datapath*filename*".jld2",
                                                      overwrite_existing = true
                                                      )


# action!
run!(simulation)


# create animation
ω_timeseries = FieldTimeSeries(datapath*filename*".jld2", "Ω")
s_timeseries = FieldTimeSeries(datapath*filename*".jld2", "S")
U_timeseries = FieldTimeSeries(datapath*filename*".jld2", "U")
PV_timeseries = FieldTimeSeries(datapath*filename*".jld2", "PV")


times = ω_timeseries.times

xω, yω, zω = nodes(ω_timeseries)
xs, ys, zs = nodes(s_timeseries)
xu, yu, zu = nodes(U_timeseries)



set_theme!(Theme(fontsize=24))

@info "Making a neat movie of vorticity and speed..."

fig = Figure(resolution=(1000,1000))

axis_kwargs = (
               xlabel="x",
               ylabel="y",
               limits=((0,Lx), (0,Ly)),
               aspect=AxisAspect(1)
               )

ax_ω = Axis(fig[2,1]; title="Vorticity", axis_kwargs...)
ax_s = Axis(fig[2,2]; title="Speed", axis_kwargs...)
ax_PV = Axis(fig[3,1]; title="zonal mean ζ+βy", (xlabel="[m/s]", ylabel="y",limits=((-2,β*Lx+2), (0,Ly)))...)
ax_U = Axis(fig[3,2]; title="zonal mean u", (xlabel="[1/s]", ylabel="y",limits=((-0.2,0.2), (0,Ly)))...)

n = Observable(1)

ω = @lift interior(ω_timeseries[$n], :, :, 1)
s = @lift interior(s_timeseries[$n], :, :, 1)
U = @lift interior(U_timeseries[$n], 1, :, 1)
PV = @lift interior(PV_timeseries[$n], 1, :, 1).+yω*β


heatmap!(ax_ω, xω, yω, ω; colormap=:balance, colorrange=(-2,2))
heatmap!(ax_s, xs, ys, s; colormap=:speed, colorrange=(0,0.2))
lines!(ax_U, U, yu)
lines!(ax_PV, PV, yω)


title = @lift "t="*string(round(times[$n],digits=2))
Label(fig[1, 1:2], title, fontsize=24, tellwidth=false)

current_figure()
fig

frames = 1:length(times)

record(fig, animationpath*filename*".mp4", frames, framerate=24) do i
    n[] = i
end