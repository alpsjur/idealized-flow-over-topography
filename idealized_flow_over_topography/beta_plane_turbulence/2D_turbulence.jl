using Oceananigans
using Statistics
using CairoMakie
using Printf

Nx = 128
Ny = 128

Lx = 2π 
Ly = 2π    

# misc parameters
β = 0          # planetary beta
f = 0          # rotation 

# simulation parameters
Δt = 0.2
stop_time = 600
save_fields_interval = 0.3

grid = RectilinearGrid(size=(Nx, Ny), 
                        x = (0, Lx),
                        y = (0, Ly),
                        topology=(Periodic, Bounded, Flat)
                        )

coriolis = BetaPlane(f₀=f, β=β)

model = NonhydrostaticModel(; grid,
                            timestepper = :RungeKutta3,
                            advection = UpwindBiasedFifthOrder(),
                            closure = ScalarDiffusivity(ν=1e-5),
                            coriolis = coriolis
                            )

u, v, w = model.velocities

uᵢ = rand(size(u)...)
vᵢ = rand(size(v)...)

uᵢ .-= mean(uᵢ)
vᵢ .-= mean(vᵢ)

# Two-dimensional Taylor Green Vortex
#m = 2    # number of vortecies
#uᵢ(x,y,z) = -cos(x*m)*sin(y*m)*0.2
#vᵢ(x,y,z) = sin(x*m)*cos(y*m)*0.2

set!(model, u=uᵢ, v=vᵢ)

simulation = Simulation(model, Δt=Δt, stop_time=stop_time)

# logging simulation progress
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 5.1f, wall time: % 15s, max |u|: % 5.3f\n",
                        sim.model.clock.iteration,
                        sim.model.clock.time,
                        prettytime(1e-9 * (time_ns() - start_time)),
                        maximum(abs, sim.model.velocities.u),
                        )
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

u, v, w = model.velocities

ω = ∂x(v) - ∂y(u)

s = sqrt(u^2+v^2)

U = Average(u, dims=1)
PV = Average(ω, dims=1)

#filename = "2D_turbulence_test"
filename = "2D_turbulence_β=$β"
datapath = "idealized_flow_over_topography/beta_plane_turbulence/data/"
animationpath = "idealized_flow_over_topography/beta_plane_turbulence/animations/"

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; ω, s, u, v, U, PV),
                                                      schedule = TimeInterval(save_fields_interval),
                                                      filename = datapath*filename*".jld2",
                                                      overwrite_existing = true
                                                      )


run!(simulation)

ω_timeseries = FieldTimeSeries(datapath*filename*".jld2", "ω")
s_timeseries = FieldTimeSeries(datapath*filename*".jld2", "s")
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
ax_PV = Axis(fig[3,1]; title="zonal mean ζ+βy", (xlabel="[m/s]", ylabel="y",limits=((-2,β*Lx+2), (0,Lx)))...)
ax_U = Axis(fig[3,2]; title="zonal mean u", (xlabel="[1/s]", ylabel="y",limits=((-0.2,0.2), (0,Lx)))...)

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