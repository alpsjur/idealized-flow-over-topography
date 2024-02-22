import plotly.graph_objects as go
import numpy as np
import xarray as xr

ds = xr.open_dataset("synthetic_data.nc")

s = 2
sizeref = 1
opacity = 0.6

x = ds.x.values[::s]
y = ds.y.values[::s]
z = ds.z.values[::s]

u = ds.u.values[::s,::s,::s,:]
v = ds.v.values[::s,::s,::s,:]
w = ds.w.values[::s,::s,::s,:]

t = ds.t.values

X, Y, Z = np.meshgrid(x, y, z)


# Number of time steps
num_steps = len(t)


time_index = 0
u_slice = u[:,:,:,time_index].flatten()
v_slice = v[:,:,:,time_index].flatten()
w_slice = w[:,:,:,time_index].flatten()


# Determine global velocity magnitude limits for consistent color scaling
velocity_magnitude = np.sqrt(u**2 + v**2 + w**2)
cmin, cmax = velocity_magnitude.min(), velocity_magnitude.max()


# Create the figure and add the first cone trace
fig = go.Figure(data=go.Cone(x=X.flatten(), y=Y.flatten(), z=Z.flatten(),
                             u=u_slice.flatten(), v=v_slice.flatten(), w=w_slice.flatten(),
                             cmin = cmin, cmax = cmax,
                             sizemode="absolute", 
                             #sizemode="scaled",
                             sizeref=sizeref, 
                             opacity=opacity)
                )

# Set the layout for the 3D scene (adjust axis limits as needed)
fig.update_layout(scene=dict(#aspectratio=dict(x=1, y=1, z=0.8),
                             xaxis=dict(range=[X.min(), X.max()]),
                             yaxis=dict(range=[Y.min(), Y.max()]),
                             zaxis=dict(range=[Z.min(), Z.max()])
                             )
                  )


frames = []

for t in range(num_steps):  # Assuming the first dimension of u is time
    frames.append(go.Frame(data=[go.Cone(x=X.flatten(), y=Y.flatten(), z=Z.flatten(),
                                         u=u[:,:,:,t].flatten(), v=v[:,:,:,t].flatten(), w=w[:,:,:,t].flatten(),
                                         sizemode="absolute", 
                                         #sizemode="scaled",
                                         sizeref=sizeref, 
                                         opacity=opacity,
                                         cmin = cmin, cmax = cmax,
                                         )],
                           name=str(t)))

fig.frames = frames

fig.update_layout(
    updatemenus=[{
        "type": "buttons",
        "buttons": [
            {
                "label": "Play",
                "method": "animate",
                "args": [None, {"frame": {"duration": 500, "redraw": True}, "fromcurrent": True}],
            },
            {
                "label": "Pause",
                "method": "animate",
                "args": [[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}],
            },
        ],
        "direction": "left",
        "pad": {"r": 10, "t": 87},
        "showactive": False,
        "x": 0.1,
        "xanchor": "right",
        "y": 0,
        "yanchor": "top"
    }]
)


# Continue from the previous setup

# Define slider steps
slider_steps = []
for t in range(len(frames)):
    step = {
        "method": "animate",
        "args": [[frames[t].name], {"frame": {"duration": 500, "redraw": True}, "mode": "immediate", "transition": {"duration": 300}}],
        "label": str(t),
    }
    slider_steps.append(step)

# Add sliders to the layout
fig.update_layout(
    sliders=[{
        "pad": {"b": 10, "t": 60},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": slider_steps,
        "currentvalue": {
            "visible": True,
            "prefix": "Time step: ",
            "xanchor": "right"
        },
        "transition": {"duration": 300},
    }]
)

# Update the layout for play/pause buttons if not already done
fig.update_layout(
    updatemenus=[{
        "type": "buttons",
        "buttons": [
            {
                "label": "Play",
                "method": "animate",
                "args": [None, {"frame": {"duration": 500, "redraw": True}, "fromcurrent": True, "transition": {"duration": 300}}],
            },
            {
                "label": "Pause",
                "method": "animate",
                "args": [[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate", "transition": {"duration": 0}}],
            },
        ],
        "direction": "left",
        "pad": {"r": 10, "t": 87},
        "showactive": False,
        "x": 0.1,
        "xanchor": "right",
        "y": 0,
        "yanchor": "top"
    }]
)

fig.write_html("cone.html")

# Show the figure
fig.show()



