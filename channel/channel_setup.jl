# Import necessary packages
using Oceananigans              # Main package for ocean simulation
using Oceananigans.Units        # For physical units (e.g., meters, seconds)
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, PartialCellBottom  # For handling complex geometries
using Random                    # For generating random numbers
using Printf                    # For formatted output
using CairoMakie                # For visualization
using CUDA

figurepath = "channel/figures/"

# Simulation parameters
Δt = 30second                   # Time step size
stop_time = 400days             # Simulation stop time
save_fields_interval = 1day     # Interval for saving output fields
average_window = 1day           # Averaging window for output data

# Grid parameters
Lx = 416kilometers              # Domain length in x-direction
Ly = 1024kilometers             # Domain length in y-direction
Lz = 2300meters                 # Domain depth
dx = 2kilometers                # Grid spacing in x-direction
#dy = 2kilometers                # Grid spacing in y-direction
dy_slope = 1kilometer
dy_flat  = 2kilometers
Nx = Int(Lx/dx)                 # Number of grid cells in x-direction
#Ny = Int(Ly/dy)                 # Number of grid cells in y-direction
Nz = 50                         # Number of grid cells in z-direction

# Bathymetry parameters (Nummelin & Isachsen, 2024)
W  = 200kilometers               # Width parameter for bathymetry
YC = 200kilometers              # Center y-coordinate for bathymetry features
DS = 2000meters                 # Depth scale
DB = 250meters                  # Base depth
σ  = 10meters                   # Standard deviation for random noise in topography

# Forcing parameters
const τ = -0.05/1000                  # Wind stress (kinematic forcing)

# Bottom friction
const Cd = 0.002                      # Quadratic drag coefficient []

# Rotation
f = 1e-4
coriolis = FPlane(f)


decay_from_LR(bmax, LR, f) = ((LR*abs(f)/(2))^2)/bmax


# Stratification parameters
bmax = 3e-2                                    # Maximum buoyancy anomaly
LR = 35kilometers                           # Deformation radius over deep ocean
decay = decay_from_LR(bmax, LR, f)          # Decay scale for bouyancy profile [m]


# Turbulence closures parameters for vertical and horizontal mixing
κh = 100    # [m²/s] horizontal diffusivity (tracers)
νh = 100    # [m²/s] horizontal viscocity   (momentum)
κz = 1e-2   # [m²/s] vertical diffusivity
νz = 1e-2   # [m²/s] vertical viscocity
vertical_closure = VerticalScalarDiffusivity(ν = νz, κ = κz)  
#vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = κz)                  
horizontal_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)
closure = (horizontal_closure, vertical_closure)

#closure = ScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=νz, κ=κz)



# Run on GPU (wow, fast!) if available. Else run on CPU
if CUDA.functional()
    architecture = GPU()
    @info "Running on GPU"
else
    architecture = CPU()
    @info "Running on CPU"
end


# Stretching function for vertical coordinates, with increased resolution near surface
function Σ(nk)     
    refinement = 0.2 # fraction of points for surface layer
    depth = 0.1      # fractional depth of surface layer
    rate = 11        # tuned manually for smooth change in Δz
    if nk < 1-refinement
        return (1-depth)*nk/(1-refinement)
    else 
        return (1-depth) + (tanh((nk-1+refinement)*rate)/tanh(refinement*rate))*depth
    end
end

# Function for generating z-coordinate faces,
z_faces(k) = Lz * (Σ(k/Nz)-1)   # Uses stretching function Σ for vertical grid spacing

# define bathymetry
function hᵢ(x, y)
    if y < (YC + W)                # shelf
        h =  -DB - 0.5*DS*(1+tanh.(π*(y-YC)/W))
    elseif Ly - y < (YC + W)       # slope
        h = -DB - 0.5*DS*(1+tanh.(π*(Ly-y-YC)/W))
    else                           # central basin
        h =  -DB - DS
    end
    # add random noise
    #h += randn()*σ
    return h
end

hᵢ(y) = hᵢ(1, y)


# Loop for creating an array of y-coordinate faces, with increased resolution over slope to resolve bathymetry
y_faces = [0]
yj = 0
while yj < Ly
    if YC-0.5*W < yj < YC+0.5*W          
        dyj = dy_slope
    elseif YC-0.5*W < Ly-yj < YC+0.5*W   
        dyj = dy_slope
    else                       
        dyj = dy_flat
    end
    global yj += dyj
    push!(y_faces, yj)
end        
Ny = length(y_faces)-1

# spesify bottom drag 
drag_u(x, y, t, u, v, Cd) = -Cd*√(u^2+v^2)*u
drag_v(x, y, t, u, v, Cd) = -Cd*√(u^2+v^2)*v

drag_u(y, t, u, v, Cd) = drag_u(1, y, t, u, v, Cd)
drag_v(y, t, u, v, Cd) = drag_v(1, y, t, u, v, Cd)

# spesify drag on immersed boundary, note different arguments from bottom drag above
immersed_drag_u(x, y, z, t, u, v, Cd) = -Cd*√(u^2+v^2)*u
immersed_drag_v(x, y, z, t, u, v, Cd) = -Cd*√(u^2+v^2)*v

immersed_drag_u(y, z, t, u, v, Cd) = immersed_drag_u(1, y, z, t, u, v, Cd)
immersed_drag_v(y, z, t, u, v, Cd) = immersed_drag_v(1, y, z, t, u, v, Cd)

# create bottom boundary conditions
drag_u_bc = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=Cd)
drag_v_bc = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=Cd)

immersed_drag_u_bc = FluxBoundaryCondition(immersed_drag_u, field_dependencies=(:u, :v), parameters=Cd)
immersed_drag_v_bc = FluxBoundaryCondition(immersed_drag_v, field_dependencies=(:u, :v), parameters=Cd)

immersed_u_bc = ImmersedBoundaryCondition(bottom = immersed_drag_u_bc)
immersed_v_bc = ImmersedBoundaryCondition(bottom = immersed_drag_v_bc)


# spesify surface forcing. Gradualy increase forcing over 10 days, then constant
τx(x, y, t, τ) = τ*tanh(t/(10days))
τy(x, y, t, τ) = 0

τx(y, t, τ) = τx(1, y, t, τ)
τy(y, t, τ) = τy(1, y, t, τ)

# create surface boundary conditions
τx_bc = FluxBoundaryCondition(τx, parameters=τ)
τy_bc = FluxBoundaryCondition(τy, parameters=τ)


# collect boundary conditions
u_bc = FieldBoundaryConditions(
                                bottom=drag_u_bc, 
                                immersed=immersed_u_bc, 
                                top=τx_bc
                                )
v_bc = FieldBoundaryConditions(
                                bottom=drag_v_bc, 
                                immersed=immersed_v_bc, 
                                top=τy_bc
                                )



# linear equation of state, for finding buoyancy profile from T
buoyancy_from_T(T) = 9.80665.*(1.67e-4.*T .- 7.80e-4.*34)

# initial buoyancy profiles  #TODO sjekk dette Finn decay scale fra Ld i ALeksis artikkel
#initial_temperature(x, y, z) = Tmax.*exp.(z./decay)
initial_buoyancy(x, y, z) = bmax.*exp.(z./decay)            #  [m/s2]
initial_buoyancy(y, z) = initial_buoyancy(1, y, z)


