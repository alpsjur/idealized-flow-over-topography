# Import required packages for visualization, ocean simulation, formatted output, and data handling
using CairoMakie
using Oceananigans
using Printf
using JLD2

# Define the path to the saved output file containing simulation data
filename = "channel/data/3D_channel_nostrat.jld2"

file = jldopen(filename)
xindex = 1

# Open the JLD2 file and extract time series data 
vv_timeseries = FieldTimeSeries(filename, "v")
ww_timeseries = FieldTimeSeries(filename, "w")
b_timeseries = FieldTimeSeries(filename, "b")
η_timeseries = FieldTimeSeries(filename, "η")


# Extract time points and bottom height 
times = b_timeseries.times
h = b_timeseries.grid.immersed_boundary.bottom_height
h = interior(h,xindex,:,1)  # Adjust the bottom height array for visualization

# Get coordinate arrays 
#xv, yv, zv = nodes(v_timeseries[1])  
xc, yc, zc = nodes(b_timeseries[1]) 

# Shift u and v to same grid point, calculate speed
v_timeseries = deepcopy(b_timeseries)
w_timeseries = deepcopy(b_timeseries)
s_timeseries = deepcopy(b_timeseries)

ystep = 20
zstep = 5

for i in 1:length(times)
    vᵢ = vv_timeseries[i]
    wᵢ = ww_timeseries[i]

    v_timeseries[i] .= @at (Center, Center, Center) vᵢ
    w_timeseries[i] .= @at (Center, Center, Center) wᵢ
    s_timeseries[i] .= @at (Center, Center, Center) sqrt(vᵢ^2+wᵢ^2)
end


# Initialize logging for the animation creation process
@info "Making an animation from saved data..."

# Create an observable integer for indexing through time series data
n = Observable(1)

# Define a title for the animation with the time variable dynamically updated
title = @lift @sprintf("%20s, i = %i", prettytime(times[$n]), xindex)

# Extract the interior data for the u and b fields at the current time step, dynamically updated
vₙ = @lift interior(v_timeseries[$n], xindex, 1:ystep:length(yc), 1:zstep:length(zc))
wₙ = @lift interior(w_timeseries[$n], xindex, 1:ystep:length(yc), 1:zstep:length(zc))
bₙ = @lift interior(b_timeseries[$n], xindex, :, :)
ηₙ = @lift interior(η_timeseries[$n], xindex, :)
sₙ = @lift interior(s_timeseries[$n], xindex, :, :)

# Set limits for the velocity color scale
slim = maximum(interior(s_timeseries))

bmin = minimum(interior(b_timeseries))
bmax = maximum(interior(b_timeseries))

ηmax = maximum(interior(η_timeseries))
ηmin = minimum(interior(η_timeseries))

# Define common axis keywords for both plots
axis_kwargs = (
               #xlabel = "Cross-channel distance [km]",
               ylabel = "Depth [m]",
               limits = ((0, 1024*1e3), (-2250, 10)),
               titlesize = 20
               )

# Create a figure object for the animation
fig = Figure(size = (1200, 1100))

# Create axes 
ax_η = Axis(fig[2, 1]; 
    title = "Sea surface height", 
    limits = ((0, 1024*1e3), ((ηmin*1.1, ηmax*1.1))),
    titlesize = 20,
    ylabel = "η [m]"
    )
ax_v = Axis(fig[3:4, 1]; 
    title = "velocity [m/s]", 
    axis_kwargs...
    )
ax_b = Axis(fig[5:6, 1]; 
    title = "buoyancy", 
    xlabel = "Cross-channel distance [m]", 
    axis_kwargs...
    )

# Add a title 
fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

# Create line plot of sea surface height
lines!(ax_η, yc, ηₙ, color=:black, linewidth=3)

# Create a heatmap for the velocity field
hm_s = heatmap!(ax_v, yc, zc, sₙ; colorrange = (0, slim), colormap = :speed)
ar_u = arrows!(ax_v, yc[1:ystep:end], zc[1:zstep:end], vₙ, wₙ, 
    lengthscale = 5e7,
    #arrowcolor=sₙ, linecolor=sₙ, colorrange = (0, slim), colormap = :speed
)
Colorbar(fig[3:4, 2], hm_s)
band!(ax_v, yc, minimum(h), h, alpha=0.5, color=:lightgray)


# Create a heatmap for the buoyancy field
hm_b = heatmap!(ax_b, yc, zc, bₙ; colorrange = (bmin, bmax), colormap = :dense)
Colorbar(fig[5:6, 2], hm_b)
band!(ax_b, yc, minimum(h), h, alpha=0.5, color=:lightgray)

# Define the frame range for the animation
frames = 1:length(times)

# Record the animation, updating the figure for each time step
CairoMakie.record(fig, "channel/animations/3D_channel_nonstrat_slice.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")  # Log progress without creating a new line for each frame
    n[] = i             # Update the observable to the current frame index
end
