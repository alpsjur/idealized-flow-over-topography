import feather as f
import xarray as xr

file = "2D_channel_nostrat"
datapath = "topographies/data/"

x = f.read_dataframe(datapath+"x.feather").values[0]
y = f.read_dataframe(datapath+"y.feather").values[:,0]
h = f.read_dataframe(datapath+"h.feather").values

ds = xr.Dataset(
    data_vars = dict(
        h =  (["x", "y"], h),
    ),
    coords = dict(
        x = x,
        y = y
    )
)

ds.to_netcdf(datapath+file+".nc")