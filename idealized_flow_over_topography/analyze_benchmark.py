import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("benchmark_results.txt", sep=" ", header=0)

# data processing

# calculate number of gridpoints
df["gridpoints"] = df["Nx"]*df["Ny"]*df["Nz"]
df["time [hours]"] = df["time"]/(60*60)


# organizing data for plotting

# include only data with same number of time steps
steps = 1000
df = df.loc[df["timesteps"]==steps]

# divide between architecture
dfCPU = df.loc[df["architecture"]=="CPU"]
dfGPU = df.loc[df["architecture"]=="GPU"]

# divide betweennumber of threads
dfGPUthread1 = dfGPU.loc[dfGPU["threads"]==1]
dfGPUthread4 = dfGPU.loc[dfGPU["threads"]==4]

dfCPUthread1 = dfCPU.loc[dfCPU["threads"]==1]
dfCPUthread4 = dfCPU.loc[dfCPU["threads"]==4]


# plot time as function of gridpoints
fig, ax = plt.subplots()

dfCPUthread1.plot.scatter(x="gridpoints", y="time [hours]", ax=ax, 
                    color="red", 
                    marker = "o",
                    label="CPU 1 tread"
                    )
dfCPUthread4.plot.scatter(x="gridpoints", y="time [hours]", ax=ax, 
                    color="red", 
                    marker = "x",
                    label="CPU 4 treads"
                    )
dfGPUthread1.plot.scatter(x="gridpoints", y="time [hours]", ax=ax, 
                    color="blue", 
                    marker = "o",
                    label="GPU 1 tread"
                    )
dfGPUthread4.plot.scatter(x="gridpoints", y="time [hours]", ax=ax, 
                    color="blue", 
                    marker = "x",
                    label="GPU 4 treads"
                    )

ax.set_xscale("log")
ax.set_title(f"timesteps = {steps}")
ax.legend()

# temporary
fig.savefig("benchmark1.png")



# pick out one horizontal resolution, see dependenze on number of z-levels
Nh = 128

dfGPUnz = dfGPU.loc[dfGPU["Nx"]==Nh]
dfCPUnz = dfCPU.loc[dfCPU["Nx"]==Nh]

fig, ax = plt.subplots()

dfCPUnz.plot.scatter(x="Nz", y="time [hours]", ax=ax, 
                    color="red", 
                    label="CPU"
                    )

dfGPUnz.plot.scatter(x="Nz", y="time [hours]", ax=ax, 
                    color="blue", 
                    label="GPU"
                    )
ax.set_title(f"Nx = Ny = {Nh}, timesteps = {steps}")
ax.legend()

# temporary
fig.savefig("benchmark2.png")

plt.show()

