import plotly.graph_objects as go
import numpy as np
import xarray as xr
import os

filename = "3D_channel_nostrat"
datapath = "topographies/data/"
ds = xr.open_dataset(datapath+filename+".nc")

x = ds.x*1e-3
y = ds.y*1e-3
t = ds.t  # Assuming 't' is your time coordinate

X, Y = np.meshgrid(x, y)
Z = ds.h.values.T  # Initial surface

# Initialize figure with the static surface
fig = go.Figure(data=[go.Surface(z=Z, x=X, y=Y, colorscale="Viridis", name="Base")])

frames = []

for time_step in t:
    Z_eta = ds.eta.sel(t=time_step).values.T  # Adjust indexing as necessary for your data
    U = ds.u.sel(t=time_step).values.T  # Adjust indexing as necessary
    V = ds.v.sel(t=time_step).values.T  # Adjust indexing as necessary
    
    # Quiver plot data preparation
    # Normally you would downsample U, V for clarity using slicing: U[::skip], V[::skip]
    skip = 5  # Skip every 5 points for clearer quiver plot
    quiver_x = X[::skip, ::skip].flatten()
    quiver_y = Y[::skip, ::skip].flatten()
    quiver_u = U[::skip, ::skip].flatten()
    quiver_v = V[::skip, ::skip].flatten()
    
    # Calculate arrow endpoints
    quiver_x_end = quiver_x + quiver_u
    quiver_y_end = quiver_y + quiver_v
    
    # Quiver plot (vector field) as a series of arrows
    quiver_traces = []
    for i in range(len(quiver_x)):
        quiver_traces.append(go.Scatter3d(x=[quiver_x[i], quiver_x_end[i]], 
                                          y=[quiver_y[i], quiver_y_end[i]], 
                                          z=[Z[0,0], Z[0,0]],  # Place arrows at the base surface level or adjust as needed
                                          mode='lines', 
                                          line=dict(color='black', width=2)))
    
    # Create a frame for the current time step
    frame = go.Frame(data=[go.Surface(z=Z_eta, x=X, y=Y, colorscale="Viridis")] + quiver_traces,
                     name=str(time_step.values))  # Use the time step value as the frame name
    frames.append(frame)

# Add frames to the figure
fig.frames = frames

# Configure animation
fig.update_layout(
    updatemenus=[dict(type="buttons",
                      showactive=False,
                      y=1,
                      x=0.8,
                      xanchor="left",
                      yanchor="bottom",
                      pad=dict(t=45, r=10),
                      buttons=[dict(label="Play",
                                    method="animate",
                                    args=[None, dict(frame=dict(duration=500, redraw=True), 
                                                     fromcurrent=True, 
                                                     mode="immediate")])])],
    sliders=[dict(steps=[dict(method='animate', args=[[f.name], 
                                                       dict(mode="immediate", frame=dict(duration=500, redraw=True))], 
                              label=f.name) for f in frames],
                   transition=dict(duration=0),
                   x=0,
                   y=0,
                   currentvalue=dict(font=dict(size=12), prefix="Time: ", visible=True),
                   len=1.0)]
)

fig.show()
