using JLD2
using Feather
using  DataFrames
using Oceananigans

datapath = "channel/data/"
filename = "2D_channel_nostrat.jld2"

outpath = "topographies/data/"

# read time series to extract grid from
ts = FieldTimeSeries(datapath*filename, "b")

# Get bottom height array
h = ts.grid.immersed_boundary.bottom_height
h = collect(interior(h,:,:,1))

# Get coordinate arrays  
x, y, z = nodes(ts[1]) 
x = collect(x)
y = collect(y)

# Create DataFrames 
dfx = DataFrame(x=x)
dfy = DataFrame(y=y)
dfh = DataFrame(h, :auto)

# Write data in Feather format 
Feather.write(outpath*"x.feather", dfx)
Feather.write(outpath*"y.feather", dfy)
Feather.write(outpath*"h.feather", dfh)