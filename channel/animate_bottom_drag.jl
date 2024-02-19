# Import required packages for visualization, ocean simulation, formatted output, and data handling
using CairoMakie
using Oceananigans
using Printf
using JLD2

# Define the path to the saved output file containing simulation data
filename = "channel/data/test_bottom_drag_output.jld2"


# Open the JLD2 file and extract time series data 
u_bc_ts = FieldTimeSeries(filename, "u_bc_op")
v_bc_ts = FieldTimeSeries(filename, "v_bc_op")
u_im_bc_ts = FieldTimeSeries(filename, "u_im_bc_op")
v_im_bc_ts = FieldTimeSeries(filename, "v_im_bc_op")

xu, yu = nodes(u_bc_ts[1])
xv, yv = nodes(v_bc_ts[1])

# Extract time points and bottom height 
times = u_bc_ts.times
h = u_bc_ts.grid.immersed_boundary.bottom_height
h = interior(h,1,:,1)  # Adjust the bottom height array for visualization


# Initialize logging for the animation creation process
@info "Making an animation from saved data..."

# Create an observable integer for indexing through time series data
n = Observable(1)

# Define a title for the animation with the time variable dynamically updated
title = @lift @sprintf("%20s", prettytime(times[$n]))

# Extract the interior data for the current time step, dynamically updated
xi = 10
u_bcₙ = @lift interior(u_bc_ts[$n], xi, :)
v_bcₙ = @lift interior(v_bc_ts[$n], xi, :)
u_im_bcₙ = @lift interior(u_im_bc_ts[$n], xi, :)
v_im_bcₙ = @lift interior(v_im_bc_ts[$n], xi, :)


# Find values to be used for axis limits
umax = maximum([maximum(interior(u_bc_ts)), maximum(interior(u_im_bc_ts))])
vmax = maximum([maximum(interior(v_bc_ts)), maximum(interior(v_im_bc_ts))])

umin = minimum([minimum(interior(u_bc_ts)), minimum(interior(u_im_bc_ts))])
vmin = minimum([minimum(interior(v_bc_ts)), minimum(interior(v_im_bc_ts))])

# Create a figure object for the animation
fig = Figure(size = (1200, 1100))

# Create axes 
ax_u = Axis(fig[2, 1]; 
    title = "u values", 
    limits = ((0, 4e4), ((umin*1.1, umax*1.1))),
    titlesize = 20,
    #ylabel = ""
    )

ax_v = Axis(fig[3, 1]; 
    title = "v values", 
    limits = ((0, 4e4), ((vmin*1.1, vmax*1.1))),
    titlesize = 20,
    #ylabel = ""
    )

ax_h = Axis(fig[4, 1]; 
    title = "bottom height", 
    limits = ((0, 4e4), ((minimum(h), 0))),
    titlesize = 20,
    xlabel = "y"
    #ylabel = ""
    )


# Add a title 
fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

# plot lines 
vlines!(ax_u, yu, color = :gray)
lines!(ax_u, yu, u_bcₙ; linewidth = 4, label="u_bc_op")
lines!(ax_u, yu, u_im_bcₙ; linewidth = 4, linestyle = :dash, label="u_im_bc_op")

vlines!(ax_v, yv, color = :gray)
lines!(ax_v, yv, v_bcₙ; linewidth = 4, label="v_bc_op")
lines!(ax_v, yv, v_im_bcₙ; linewidth = 4, linestyle = :dash, label="v_im_bc_op")

vlines!(ax_h, yu, color = :gray)
band!(ax_h, yu, minimum(h), h, alpha=0.5, color=:lightgray)

# Define the frame range for the animation
frames = 1:length(times)

# Record the animation, updating the figure for each time step
CairoMakie.record(fig, "channel/animations/test_bottom_drag_output.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")  # Log progress without creating a new line for each frame
    n[] = i             # Update the observable to the current frame index
end