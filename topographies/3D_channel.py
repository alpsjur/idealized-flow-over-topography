import plotly.graph_objects as go
import numpy as np
import xarray as xr 

datapath = "channel/data/"
filename = "2D_channel_nostrat.jld2"

ds = xr.open_dataset(datapath+filename, engine='h5netcdf')