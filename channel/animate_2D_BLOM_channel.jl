# Import required packages for visualization, ocean simulation, formatted output, and data handling
using CairoMakie
using Oceananigans
using Printf
using JLD2

# Define the path to the saved output file containing simulation data
filename = "channel/data/2D_BLOM_channel.jld2"

file = jldopen(filename)

# Open the JLD2 file and extract time series data 
u_timeseries = FieldTimeSeries(filename, "u")
b_timeseries = FieldTimeSeries(filename, "b")

# Extract time points and bottom height 
times = b_timeseries.times
h = b_timeseries.grid.immersed_boundary.bottom_height
h = interior(h,1,:,1)  # Adjust the bottom height array for visualization

# Get coordinate arrays 
xu, yu, zu = nodes(u_timeseries[1])  
xc, yc, zc = nodes(b_timeseries[1])  

# Initialize logging for the animation creation process
@info "Making an animation from saved data..."

# Create an observable integer for indexing through time series data
n = Observable(1)

# Define a title for the animation with the time variable dynamically updated
title = @lift @sprintf("%20s", prettytime(times[$n]))

# Extract the interior data for the u and b fields at the current time step, dynamically updated
uₙ = @lift interior(u_timeseries[$n], 1, :, :)
bₙ = @lift interior(b_timeseries[$n], 1, :, :)

# Set limits for the velocity color scale
ulim = maximum(abs, interior(u_timeseries))

# Define common axis keywords for both plots
axis_kwargs = (xlabel = "Cross-channel distance [km]",
               ylabel = "Depth [m]",
               limits = ((0, 1024), (-2250, 0)),
               titlesize = 20)

# Create a figure object for the animation
fig = Figure(size = (1200, 1100))

# Create axes for velocity and buoyancy 
ax_u = Axis(fig[2, 1]; title = "u velocity [m/s]", axis_kwargs...)
ax_b = Axis(fig[3, 1]; title = "buoyancy", axis_kwargs...)

# Add a title 
fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

# Create a heatmap for the velocity field
hm_u = heatmap!(ax_u, yu*1e-3, zu, uₙ; colorrange = (-ulim, ulim), colormap = :balance)
Colorbar(fig[2, 2], hm_u)
band!(ax_u, yc*1e-3, minimum(h), h, alpha=0.5, color=:lightgray)

# Create a heatmap for the buoyancy field
hm_b = heatmap!(ax_b, yc*1e-3, zc, bₙ; colorrange = (-1, 0), colormap = :dense)
Colorbar(fig[3, 2], hm_b)
band!(ax_b, yc*1e-3, minimum(h), h, alpha=0.5, color=:lightgray)

# Define the frame range for the animation
frames = 1:length(times)

# Record the animation, updating the figure for each time step
record(fig, "channel/animations/2D_BLOM_channel.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")  # Log progress without creating a new line for each frame
    n[] = i             # Update the observable to the current frame index
end
