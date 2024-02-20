import plotly.graph_objects as go
import numpy as np
import xarray as xr
import os

filename = "3D_channel_nostrat"

os.system("julia --project=.. grid_jld2_to_feather.jl "+filename)
os.system("python grid_feather_to_nc.py "+filename)

datapath = "topographies/data/"
ds = xr.open_dataset(datapath+filename+".nc")

x = ds.x*1e-3
y = ds.y*1e-3

X, Y = np.meshgrid(x, y)
Z = ds.h.values.T

fig = go.Figure(data=[go.Surface(x=X, y=Y,z=Z,
                                 colorscale="deep_r"
                                 )])
fig.update_traces(contours_z=dict(show=True, 
                                  usecolormap=True, 
                                  #highlightcolor="limegreen", 
                                  project_z=True)
                  )


fig.update_layout(title=filename, 
                  #template="plotly_dark",
                  template="plotly",
                  scene=dict(
        xaxis=dict(range=[min(x), max(x)], title='x [km]'),  
        yaxis=dict(range=[min(y), max(y)], title='y [km]'),  
        zaxis=dict(range=[min(Z.flatten())*1.1, 0], title='Depth [m]') 
    )
)

fig.show()