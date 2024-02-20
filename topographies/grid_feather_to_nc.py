import feather as f
import xarray as xr
import sys

file = sys.argv[1]
#file = "3D_channel_nostrat"
datapath = "data/"
#datapath = "topographies/data/"

x = f.read_dataframe(datapath+"x.feather").values.flatten()
y = f.read_dataframe(datapath+"y.feather").values.flatten()
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