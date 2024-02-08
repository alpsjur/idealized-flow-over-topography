# Import necessary packages
using Oceananigans              # Main package for ocean simulation
using Oceananigans.Units        # For physical units (e.g., meters, seconds)
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, PartialCellBottom  # For handling complex geometries
using Random                    # For generating random numbers
using Printf                    # For formatted output
using CairoMakie                # For visualization

# Simulation parameters
Δt = 30second                   # Time step size
stop_time = 10days              # Simulation stop time
save_fields_interval = 1hour    # Interval for saving output fields
average_window = 1hour          # Averaging window for output data

# Grid parameters
Lx = 416kilometers              # Domain length in x-direction
Ly = 1024kilometers             # Domain length in y-direction
Lz = 2250meters                 # Domain depth
dx = 2kilometers                # Grid spacing in x-direction
dy = 2kilometers                # Grid spacing in y-direction
Nx = Int(Lx/dx)                 # Number of grid cells in x-direction
Ny = Int(Ly/dy)                 # Number of grid cells in y-direction
Nz = 51                         # Number of grid cells in z-direction

# Function for generating z-coordinate faces
z_faces(k) = Lz * (Σ(k/Nz)-1)   # Uses stretching function Σ for vertical grid spacing

# Bathymetry parameters (Nummelin & Isachsen, 2024)
W  = 75kilometers               # Width parameter for bathymetry
YC = 150kilometers              # Center y-coordinate for bathymetry features
DS = 2000meters                 # Depth scale
DB = 250meters                  # Base depth
σ  = 10meters                   # Standard deviation for random noise in topography

# Forcing parameters
τ = -0.05/1000                  # Wind stress (kinematic forcing)

# Bottom friction
R = 0.002                       # Quadratic drag coefficient

# Stratification parameters
Tmax = 15                       # Maximum temperature
bmax = 0.76-1                   # Maximum buoyancy anomaly
decay = 100meter                # Decay scale for temperature and buoyancy profiles


# Turbulence closures parameters for vertical and horizontal mixing
κh = 0.5e-5 # [m²/s] horizontal diffusivity
νh = 30.0   # [m²/s] horizontal viscocity
κz = 0.5e-5 # [m²/s] vertical diffusivity
νz = 3e-4   # [m²/s] vertical viscocity
vertical_closure = ScalarDiffusivity(ν = νz, κ = κz)                  
horizontal_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)


# Rotation
coriolis = FPlane(1e-4)


# Stretching function for vertical coordinates
function Σ(nk)     
    refinement = 0.2 # fraction of points for surface layer
    depth = 0.1      # fractional depth of surface layer
    rate = 12        # tuned manually for smooth change in Δz
    if nk < 1-refinement
        return (1-depth)*nk/(1-refinement)
    else 
        return (1-depth) + (atan((nk-1+refinement)*rate)/atan(refinement*rate))*depth
    end
end



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

# spesify bottom drag 
drag_u(x, y, t, u, v, R) = -R*u*√(u^2*v^2)
drag_v(x, y, t, u, v, R) = -R*v*√(u^2*v^2) 

drag_u(y, t, u, v, R) = drag_u(1, y, t, u, v, R)
drag_v(y, t, u, v, R) = drag_v(1, y, t, u, v, R)

# spesify drag on immersed boundary, note different arguments from bottom drag above
immersed_drag_u(x, y, z, t, u, v, R) = -R*u*√(u^2*v^2)
immersed_drag_v(x, y, z, t, u, v, R) = -R*v*√(u^2*v^2)

immersed_drag_u(y, z, t, u, v, R) = immersed_drag_u(1, y, z, t, u, v, R)
immersed_drag_v(y, z, t, u, v, R) = immersed_drag_v(1, y, z, t, u, v, R)

# create bottom boundary conditions
drag_u_bc = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=R)
drag_v_bc = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=R)

immersed_drag_u_bc = FluxBoundaryCondition(immersed_drag_u, field_dependencies=(:u, :v), parameters=R)
immersed_drag_v_bc = FluxBoundaryCondition(immersed_drag_v, field_dependencies=(:u, :v), parameters=R)

immersed_u_bc = ImmersedBoundaryCondition(bottom = immersed_drag_u_bc)
immersed_v_bc = ImmersedBoundaryCondition(bottom = immersed_drag_v_bc)


# spesify surface forcing
τx(x, y, t, τ) = τ
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

# initial buoyancy profiles
#initial_temperature(x, y, z) = Tmax.*exp.(z./decay)
initial_buoyancy(x, y, z) = (bmax+1).*exp.(z./decay).-1

initial_buoyancy(y, z) = initial_buoyancy(1, y, z)


