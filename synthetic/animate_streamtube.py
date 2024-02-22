import plotly.graph_objects as go
import numpy as np
import xarray as xr

ds = xr.open_dataset("synthetic_data.nc")

x = ds.x.values
y = ds.y.values
z = ds.z.values

u = ds.u.values
v = ds.v.values
w = ds.w.values

t = ds.t.values

X, Y, Z = np.meshgrid(x, y, z)


# Number of time steps
num_steps = len(t)

# Initial figure setup with fixed axis limits and initial layout
fig = go.Figure(layout=dict(scene=dict(
    xaxis=dict(range=[X.min(), X.max()]),  # Fixing X-axis limits
    yaxis=dict(range=[Y.min(), Y.max()]),  # Fixing Y-axis limits
    zaxis=dict(range=[Z.min(), Z.max()]),  # Fixing Z-axis limits
)))

# Determine global velocity magnitude limits for consistent color scaling
velocity_magnitude = np.sqrt(u**2 + v**2 + w**2)
cmin, cmax = velocity_magnitude.min(), velocity_magnitude.max()

# Assuming X, Y, Z are 3D arrays and are constant over time
# Flatten the arrays for streamtube plotting
X_flat = X.flatten()
Y_flat = Y.flatten()
Z_flat = Z.flatten()

s = 5
starts = dict(x=X[::s, 0, ::s].flatten(), y=Y[::s, 0, ::s].flatten(), z=Z[::s, 0, ::s].flatten())


# Loop through each time step to create a frame
for time_index in range(num_steps):
    u_slice = u[:,:,:,time_index].flatten()
    v_slice = v[:,:,:,time_index].flatten()
    w_slice = w[:,:,:,time_index].flatten()

    
    
    # Create a streamtube trace for this time step
    frame_trace = go.Streamtube(x=X_flat, y=Y_flat, z=Z_flat, 
                                u=u_slice, v=v_slice, w=w_slice, 
                                starts=starts, 
                                cmin=cmin, cmax=cmax, 
                                colorscale='Viridis', 
                                showscale=False, 
                                sizeref = 0.3,
                                )

    # Add this trace to a new frame in the figure
    fig.add_trace(frame_trace)

# Make only the first trace visible
for i in range(num_steps):
    fig.data[i].visible = False
fig.data[0].visible = True

# Create slider steps for animation control
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": f"Time step: {i}"}],
    )
    step["args"][0]["visible"][i] = True  # Toggle visibility
    steps.append(step)

# Update figure layout to include slider and fix scene aspect ratio
fig.update_layout(
    sliders=[dict(active=0, currentvalue={"prefix": "Time step: "}, pad={"t": 50}, steps=steps)],
    #updatemenus=[dict(type="buttons",
    #                      buttons=[dict(label="Play",
    #                                    method="animate",
    #                                    args=[None])])],
    title_text="Streamtube Animation Over Time",
    scene=dict(aspectmode='cube')  # Optional: Fix the aspect ratio
)

fig.write_html("streamtube.html")
fig.show()