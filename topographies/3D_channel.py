import plotly.graph_objects as go
import numpy as np
import xarray as xr 

datapath = "../channel/data"
filename = "3D_channel.jld2"

ds = xr.open_dataset(datapath+filename)