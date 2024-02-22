
"""
This discussion should be useful too: https://github.com/CliMA/Oceananigans.jl/discussions/3363
They're using the advective terms as an example, but basically it seems you can 
define all the operations you need, like du/dt, (Av*u_z)_z, and pass them to your output writer.
This way you should be able to save the momentum terms and close the budget exactly 
at any output frequency you like. You might just have to define the time derivative separately. 
That definition would have to match your simulation's time-stepping scheme, like what 
they talk about regarding the advection scheme for the advective terms.

a
"""
# time derivative of velocity 
function ∂u∂t(model) 
    # TODO fix first timestep
    Gⁿ = model.timestepper.Gⁿ.u 
    G⁻ = model.timestepper.G⁻.u 
    χ  = model.timestepper.χ.u 

    return (3/2 + χ) * Gⁿ - (1/2 + χ) * G⁻
end

function ∂v∂t(model) 
    # TODO fix first timestep
    Gⁿ = model.timestepper.Gⁿ.v 
    G⁻ = model.timestepper.G⁻.v 
    χ  = model.timestepper.χ.v 

    return (3/2 + χ) * Gⁿ - (1/2 + χ) * G⁻
end

function init_save_some_metadata!(file, model)
    file["author"] = "Anna Lina Sjur"
    return nothing
end


u, v, w = model.velocities
p = model.pressure.pHY′      # see here: https://github.com/CliMA/Oceananigans.jl/discussions/3157
η′ = model.free_surface.η
b = model.tracers.b
h = model.grid.immersed_boundary.bottom_height

η = Average(η′, dims=3)


uu = u*u 
vv = v*v
uv = u*v

ub = u*b
vb = v*b 
wb = w*b


# Info for NetCDFOutputWriter 
dims = Dict("h" => ("xC", "yC"))

output_attributes = Dict(
    "h"  => Dict("long_name" => "Bottom height", "units" => "m"),
    "eta"  => Dict("long_name" => "Free surface height", "units" => "m"),
    "p"  => Dict("long_name" => "Pressure", "units" => "N/m²"),
    "uu"  => Dict("long_name" => "Velocity advection u*u", "units" => "m²/s²"),
    "vv"  => Dict("long_name" => "Velocity advection v*v", "units" => "m²/s²"),
    "uv"  => Dict("long_name" => "Velocity advection u*v", "units" => "m²/s²"),
    "ub"  => Dict("long_name" => "Buoyancy advection u*b", "units" => "m²/s²"),
    "vb"  => Dict("long_name" => "Buoyancy advection v*b", "units" => "m²/s³"),
    "wb"  => Dict("long_name" => "Buoyancy advection w*b", "units" => "m²/s³"),
)

global_attributes = Dict(
    "author" => "Anna Lina Sjur"
)

#bottom drag
#Bottom drag, kun på bunn: https://github.com/CliMA/Oceananigans.jl/discussions/3081
#Finne verdi bare ved immersed bahtymetry: https://github.com/CliMA/Oceananigans.jl/discussions/3032

# code based on https://github.com/CliMA/Oceananigans.jl/discussions/3032 for extracting fields at immersed boundary
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node
using Oceananigans.Fields: condition_operand

@inline bottom_condition(i, j, k, grid) = immersed_peripheral_node(i, j, k-1, grid, Center(), Center(), Center())

# bottom velocities
u_b = sum(condition_operand(identity, u, bottom_condition, 0.0), dims = 3)
v_b = sum(condition_operand(identity, v, bottom_condition, 0.0), dims = 3)


# bottom pressure 
p_b = sum(condition_operand(identity, p, bottom_condition, 0.0), dims = 3)


### Code based on https://github.com/CliMA/Oceananigans.jl/discussions/3081 for extracting bottom drag
using Oceananigans.BoundaryConditions: getbc
using Oceananigans: fields

# Boundary condition extractor in "kernel function form"
@inline kernel_getbc(i, j, k, grid, boundary_condition, clock, fields) =
    getbc(boundary_condition, i, j, grid, clock, fields)

# Kernel arguments
clock = model.clock
model_fields = merge(fields(model), model.auxiliary_fields)

u_bc = u.boundary_conditions.bottom
v_bc = v.boundary_conditions.bottom

u_im_bc = u.boundary_conditions.immersed.bottom
v_im_bc = v.boundary_conditions.immersed.bottom

# Build operations. these can be passed to OutputWriters
u_bc_op = KernelFunctionOperation{Face, Center, Nothing}(kernel_getbc, grid, u_bc, clock, model_fields)
v_bc_op = KernelFunctionOperation{Center, Face, Nothing}(kernel_getbc, grid, v_bc, clock, model_fields)

u_im_bc_op = KernelFunctionOperation{Face, Center, Nothing}(kernel_getbc, grid, u_im_bc, clock, model_fields)
v_im_bc_op = KernelFunctionOperation{Face, Center, Nothing}(kernel_getbc, grid, v_im_bc, clock, model_fields)


# Build Fields
u_bc_field = Field(u_bc_op)
v_bc_field = Field(v_bc_op)

u_im_bc_field = Field(u_im_bc_op)
v_im_bc_field = Field(v_im_bc_op)


# logging simulation progress
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 15s, wall time: % 15s, max |u|: % 5.3f, max |v|: % 5.3f, max |w|: % 5.3f, max |η|: % 5.3f, next Δt: %s\n",
        sim.model.clock.iteration,
        prettytime(sim.model.clock.time),
        #sim.model.clock.time,
        prettytime(1e-9 * (time_ns() - start_time)),
        maximum(abs, u),
        maximum(abs, v),
        maximum(abs, w),
        maximum(abs, η′),
        prettytime(sim.Δt),
)
