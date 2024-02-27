# Import required packages for visualization, ocean simulation, formatted output, and data handling
using CairoMakie
using Oceananigans
using Printf
using JLD2

# Define the path to the saved output file containing simulation data
filepath = "channel/data/"
filename = "3D_channel_strat_average"

file = jldopen(filepath*filename*".jld2")


# Open the JLD2 file and extract time series data 
v_timeseries = FieldTimeSeries(filepath*filename*".jld2", "V")
w_timeseries = FieldTimeSeries(filepath*filename*".jld2", "W")
b_timeseries = FieldTimeSeries(filepath*filename*".jld2", "B")
η_timeseries = FieldTimeSeries(filepath*filename*".jld2", "H")


# Extract time points and bottom height 
times = v_timeseries.times
h = v_timeseries.grid.immersed_boundary.bottom_height
h = interior(h,1,:,1)  # Adjust the bottom height array for visualization

# Get coordinate arrays 
#xv, yv, zv = nodes(v_timeseries[1])  
xc, yc, zc = nodes(b_timeseries[1]) 

# Shift u and v to same grid point, calculate speed
vc_timeseries = deepcopy(b_timeseries)
wc_timeseries = deepcopy(b_timeseries)
s_timeseries = deepcopy(b_timeseries)

ystep = 20
zstep = 5

for i in 1:length(times)
    vᵢ = v_timeseries[i]
    wᵢ = w_timeseries[i]

    vc_timeseries[i] .= @at (Center, Center, Center) vᵢ
    wc_timeseries[i] .= @at (Center, Center, Center) wᵢ
    s_timeseries[i] .= @at (Center, Center, Center) sqrt(vᵢ^2+wᵢ^2)
end


# Initialize logging for the animation creation process
@info "Making an animation from saved data..."

# Create an observable integer for indexing through time series data
n = Observable(1)

# Define a title for the animation with the time variable dynamically updated
title = @lift @sprintf("%20s", prettytime(times[$n]))

# Extract the interior data for the u and b fields at the current time step, dynamically updated
vₙ = @lift interior(vc_timeseries[$n], 1, 1:ystep:length(yc), :)
wₙ = @lift interior(wc_timeseries[$n], 1, 1:ystep:length(yc), :)
vsₙ = @lift interior(vc_timeseries[$n], 1, 1:ystep:length(yc), 1:zstep:length(zc))
wsₙ = @lift interior(wc_timeseries[$n], 1, 1:ystep:length(yc), 1:zstep:length(zc))
bₙ = @lift interior(b_timeseries[$n], 1, :, :)
ηₙ = @lift interior(η_timeseries[$n], 1, :)
sₙ = @lift interior(s_timeseries[$n], 1, :, :)

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
ax_v = Axis(fig[3, 1]; 
    title = "velocity [m/s]", 
    ylabel = "Depth [m]",
    limits = ((0, 1024*1e3), (-25, 10)),
    titlesize = 20
    )
ax_vs = Axis(fig[4:5, 1]; 
    title = "velocity [m/s]", 
    axis_kwargs...
    )

ax_b = Axis(fig[6:7, 1]; 
    title = "buoyancy", 
    xlabel = "Cross-channel distance [m]", 
    axis_kwargs...
    )

# Add a title 
fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

# Create line plot of sea surface height
lines!(ax_η, yc, ηₙ, color=:black, linewidth=3)

# Create a heatmap for the surface velocity field
hm_s = heatmap!(ax_v, yc, zc, sₙ; colorrange = (0, slim), colormap = :speed)
ar_u = arrows!(ax_v, yc[1:ystep:end], zc, vₙ, wₙ, 
    lengthscale = 1e6,
    #arrowcolor=sₙ, linecolor=sₙ, colorrange = (0, slim), colormap = :speed
)
band!(ax_v, yc, minimum(h), h, alpha=0.5, color=:lightgray)
Colorbar(fig[3, 2], hm_s)


# Create a heatmap for the sparce velocity field
hm_s = heatmap!(ax_vs, yc, zc, sₙ; colorrange = (0, slim*0.1), colormap = :speed)
ar_u = arrows!(ax_vs, yc[1:ystep:end], zc[1:zstep:end], vsₙ, wsₙ, 
    lengthscale = 5e7,
    #arrowcolor=sₙ, linecolor=sₙ, colorrange = (0, slim), colormap = :speed
)
band!(ax_vs, yc, minimum(h), h, alpha=0.5, color=:lightgray)
Colorbar(fig[4:5, 2], hm_s)


# Create a heatmap for the buoyancy field
hm_b = heatmap!(ax_b, yc, zc, bₙ; colorrange = (bmin, bmax), colormap = :dense)
Colorbar(fig[6:7, 2], hm_b)
band!(ax_b, yc, minimum(h), h, alpha=0.5, color=:lightgray)


# Define the frame range for the animation
frames = 1:length(times)

# Record the animation, updating the figure for each time step
CairoMakie.record(fig, "channel/animations/"*filename*".mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")  # Log progress without creating a new line for each frame
    n[] = i             # Update the observable to the current frame index
end
