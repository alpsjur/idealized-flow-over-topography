using Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, PartialCellBottom
using Random
using Printf
using CairoMakie


# simulation parameters 
Δt = 120second
stop_time = 30days
save_fields_interval = 1hour
average_window = 1hour

# grid parameters 
Lx =  416kilometers
Ly = 1024kilometers
Lz = 2250meters

dx = 2kilometers
dy = 2kilometers

Nx = Int(Lx/dx)
Ny = Int(Ly/dy)
Nz = 51

# Generating function
z_faces(k) = Lz * (Σ(k/Nz)-1)


# bathymetry parameters (see Nummelin & Isachsen, 2024)
W  =   75kilometers
YC =  150kilometers
DS = 2000meters
DB =  250meters
σ  =   10meters              # std of random noise on topography


# forcing parameters
τ = -0.05/1000               # N m-2 / kg m-2  kinematic forcing

# bottom friction
R = 0.002                    # Bottom quadratic drag coeficient TODO correct number

# stratification parameters 
# svak restoring til denne profilen pga eksplisitt horisontalmiksing + numerisk diffusjon
Tmax = 15
bmax = 0.76-1
decay = 100meter 


# stretching function for z-koordinates
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




# define 3D bathymetry
function hᵢ(x, y)
    if y < (YC + W)                # shelf
        h =  -DB - 0.5*DS*(1+tanh.(π*(y-YC)/W))
    elseif Ly - y < (YC + W)       # slope
        h = -DB - 0.5*DS*(1+tanh.(π*(Ly-y-YC)/W))
    else                           # central basin
        h =  -DB - DS
    end
    # add random noise
    h += randn()*σ
    return h
end


# define 2D bathymetry
function hᵢ(y)
    if y < (YC + W)                # shelf
        h =  -DB - 0.5*DS*(1+tanh.(π*(y-YC)/W))
    elseif Ly - y < (YC + W)       # slope
        h = -DB - 0.5*DS*(1+tanh.(π*(Ly-y-YC)/W))
    else                           # central basin
        h =  -DB - DS
    end
    
    return h
end


# Turbulence closures 
κh = 0
νh = 0   
κz = 1e-5
νz = 1e-5
vertical_closure = ScalarDiffusivity(ν = νz, κ = κz)                 
horizontal_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)

# Rotation
coriolis = FPlane(1e-4)

# spesify 3D bottom drag 
drag_u(x, y, t, u, v, R) = -R*u*√(u^2*v^2)
drag_v(x, y, t, u, v, R) = -R*v*√(u^2*v^2) 

# spesify 2D bottom drag
drag_u(y, t, u, v, R) = drag_u(1, y, t, u, v, R)
drag_v(y, t, u, v, R) = drag_v(1, y, t, u, v, R)

# spesify 3D drag on immersed boundary, note different arguments from bottom drag above
immersed_drag_u(x, y, z, t, u, v, R) = -R*u*√(u^2*v^2)
immersed_drag_v(x, y, z, t, u, v, R) = -R*v*√(u^2*v^2)

# spesify 2D immersed drag
immersed_drag_u(y, z, t, u, v, R) = immersed_drag_u(1, y, z, t, u, v, R)
immersed_drag_v(y, z, t, u, v, R) = immersed_drag_v(1, y, z, t, u, v, R)

# spesify 3D surface forcing
τx(x, y, t, τ) = τ
τy(x, y, t, τ) = 0

# spesify 2D surface forcing
τx(y, t, τ) = τx(1, y, t, τ)
τy(y, t, τ) = τy(1, y, t, τ)

τx_bc = FluxBoundaryCondition(τx, parameters=τ)
τy_bc = FluxBoundaryCondition(τy, parameters=τ)


# linear equation of state, for finding buoyancy profile from T
buoyancy_from_T(T) = 9.80665.*(1.67e-4.*T .- 7.80e-4.*34)

# TODO hvilken equation of state er i Aleksis simuleringer? 
# define initial temperature TODO update

# 3D initial profiles
initial_temperature(x, y, z) = Tmax.*exp.(z./decay)
initial_buoyancy(x, y, z) = (bmax+1).*exp.(z./decay).-1

# 2D initial profile
initial_buoyancy(y, z) = initial_buoyancy(1, y, z)


