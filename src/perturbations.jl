"""
    simulate_pulse(c::Community, x, tspan)

Simulate the recovery of the community `c`
after the pulse perturbation `x` for the time span `tspan`.
In other words, it simulate the community dynamics
of initial conditions ``N^* + x`` where ``N^*`` is the vector
of species equilibrium abundances.

## Example

```julia
using Distributions
S = 5 # Number of species.
c = rand(Community, S; A_ij=Normal(0, 0.1))
x = rand(Normal(-2, 0.1), S)
u = simulate_pulse(c, x, (0, 100))
```

See also [`solve`](@ref).
"""
function simulate_pulse(c::Community, x, tspan)
    N_eq = abundance(c)
    solve(c, N_eq + x, tspan)
end
export simulate_pulse

"""
    simulate_extinctions(c::Community, idx, tspan)

Simulate the dynamics of the community `c` after the extinction
of the species of indices `idx` for the time span `tspan`.

## Example

```julia
using Distributions
c = rand(Community, 5; A_ij=Normal(0, 0.1))
simulate_extinctions(c, [1, 3], (0, 100)) # Species 1 and 3 go extinct.
```

See also [`solve`](@ref), [`simulate_pulse`](@ref).
"""
function simulate_extinctions(c::Community, idx, tspan)
    N0 = abundance(c) # Start with species equilibrium abundances...
    N0[idx] .= 0 # ... except for the extinct species.
    solve(c, N0, tspan)
end
export simulate_extinctions

"""
    simulate_press(c::Community, K_new, tspan)

Simulate the dynamics of the community `c` following a press perturbation
The press perturbation is modeled by a change in the carrying capacities.
The species carrying capacities after the perturbation are given by `K_new`.

## Example

```julia
using Distributions
c = rand(Community, 5; A_ij=Normal(0, 0.1))
K_new = c.K .- [0.9, 0, 0, 0, 0] # Lower the carrying capacity of the first species.
simulate_press(c, K_new, (0, 100))
```

See also [`solve`](@ref), [`simulate_pulse`](@ref), [`simulate_extinctions`](@ref).
"""
function simulate_press(c::Community, K_new, tspan)
    N0 = abundance(c)
    c.K .= K_new
    solve(c, N0, tspan)
end
export simulate_press

"""
    simulate_noise(c::Community, noise!::Function, tspan)

Simulate the dynamics of the community `c` with stochastic noise, around its equilibrium.
The species equilibrium abundances are given by [`abundance`](@ref).

See also [`solve`](@ref), [`simulate_pulse`](@ref), [`simulate_extinctions`](@ref), [`simulate_press`](@ref).
"""
function simulate_noise(c::Community, noise!::Function, tspan)
    N0 = abundance(c)
    solve(c, N0, tspan, noise!)
end
export simulate_noise
