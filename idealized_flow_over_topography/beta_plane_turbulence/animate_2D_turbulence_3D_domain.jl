using Oceananigans
using CairoMakie

filename = "2Dturbulence_3Ddomain"
datapath = "idealized_flow_over_topography/beta_plane_turbulence/data/"
animationpath = "idealized_flow_over_topography/beta_plane_turbulence/animations/"

# create animation
ω_timeseries = FieldTimeSeries(datapath*filename*".jld2", "Ω")
s_timeseries = FieldTimeSeries(datapath*filename*".jld2", "S")
U_timeseries = FieldTimeSeries(datapath*filename*".jld2", "U")
PV_timeseries = FieldTimeSeries(datapath*filename*".jld2", "PV")


times = ω_timeseries.times

xω, yω, zω = nodes(ω_timeseries)
xs, ys, zs = nodes(s_timeseries)


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