
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

#bottom friction force
"""
regne ut ved hjelp av likning for drag definert tidligere? Hvordan finne u ved bunn?
"""






function init_save_some_metadata!(file, model)
    file["author"] = "Anna Lina Sjur"
    return nothing
end