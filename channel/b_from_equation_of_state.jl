using Oceananigans

buoyancy = SeawaterBuoyancy()

α = buoyancy.equation_of_state.thermal_expansion
β = buoyancy.equation_of_state.haline_contraction
g = buoyancy.gravitational_acceleration

b(T, S) = g*(α*T-β*S)

T0 = 6
T1 = 0
S = 32



b0 = b(T0, S)
b1 = b(T1, S)

Δb = b0-b1