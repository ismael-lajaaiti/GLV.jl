"""
    simulate_pulse(c::Community, x, tspan)

Simulate the recovery of the community `c`
after the pulse perturbation `x` for the time span `tspan`.
In other words, it simulate the community dynamics
of initial conditions ``N^* + x`` where `N^*` is the vector
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

# Example

```julia
c = rand(Community, 5; A_ij=Normal(0, 0.1))
simulate_extinctions(c, [1, 3], (0, 100)) # Species 1 and 3 go extinct.
```

See also [`solve`](@ref), [`simulate_pulse`](@ref).
"""
function simulate_extinctions(c::Community, idx, tspan)
    N0 = abundance(c) # Start with species equilibrium abundances...
    N0[idx] .= 0 # ... expect for the extinct species.
    solve(c, N0, tspan)
end
export simulate_extinctions
