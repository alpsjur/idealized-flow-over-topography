using CairoMakie
using Oceananigans
using Oceananigans.Fields
using Oceananigans.AbstractOperations: volume

# based on https://clima.github.io/OceananigansDocumentation/stable/generated/horizontal_convection/

saved_output_filename = "data/2D_BLOM_channel_test.jld2"

# Open the file with our data
u_timeseries = FieldTimeSeries(saved_output_filename, "u")
b_timeseries = FieldTimeSeries(saved_output_filename, "b")

times = b_timeseries.times

# Coordinate arrays
xu, yu, zu = nodes(u_timeseries[1])
xc, yc, zc = nodes(b_timeseries[1])


@info "Making an animation from saved data..."

n = Observable(1)

title = @lift @sprintf("t=%1.2f", times[$n])

uₙ = @lift interior(ζ_timeseries[$n], :, 1, :)
bₙ = @lift interior(b_timeseries[$n], :, 1, :)

ulim = 0.6
blim = 0.6

axis_kwargs = (xlabel = L"Cross-channel distance [km]",
               ylabel = L"Depth [m]",
               limits = ((0, 1024), (-2250, 0)),
               titlesize = 20)

fig = Figure(size = (600, 1100))

ax_u = Axis(fig[2, 1];
            title = L"u velocity [m/s]", axis_kwargs...)

ax_b = Axis(fig[3, 1];
            title = L"buoyancy", axis_kwargs...)

fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

hm_u = heatmap!(ax_u, xu, zu, uₙ;
                colorrange = (-ulim, ulim),
                colormap = :speed)
Colorbar(fig[2, 2], hm_u)

hm_b = heatmap!(ax_b, xc, zc, bₙ;
                colorrange = (-blim, blim),
                colormap = :thermal)
Colorbar(fig[3, 2], hm_b)


frames = 1:length(times)

record(fig, "animations/2D_BLOM_channel.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end