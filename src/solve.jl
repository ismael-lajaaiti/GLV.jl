"""
    DifferentialEquations.solve(c::Community, u0, tspan; kwargs...)

Run the GLV model for community `c` with initial conditions `u0` and time span `tspan`.
The GLV model writes
```math
\\frac{\\mathrm{d} N_i}{\\mathrm{d}t} = r_i N_i \\left(\\frac{\\sum_{j\\neq i} A_{ij} N_j - N_i}{K_i}\\right)
```
where ``r`` is the growth rate, ``A`` is the interaction matrix, and ``K`` is the carrying capacity.

## Example

Two non-interacting species with different carrying capacities.

```julia
A = [-1 0; 0 -1] # Only self-interactions.
r = [1.0, 1.0]
K = [1.0, 2.0]
c = Community(A, r, K)
u0, tspan = [1.0, 1.0], (0, 10_000) # Simulation parameters.
sol = solve(c, u0, tspan) # Simulate the dynamics.
```

See also [`Community`](@ref).
"""
function solve(c::Community, u0, tspan; kwargs...)
    function f!(du, u, _, _)
        for i in eachindex(u)
            u[i] < 0 && (u[i] = 0) # Species cannot have negative abundances.
            du[i] = c.r[i] * u[i] * (1 + sum(c.A[i, :] .* u) / c.K[i])
        end
    end
    prob = ODEProblem(f!, u0, tspan)
    DifferentialEquations.solve(prob; kwargs...)
end


"""
    solve(c::Community, u0, tspan, noise!::Function; kwargs...)

Run [`solve`](@ref) with stochastic noise, given by the function `noise`.
The `noise!` function should be defined as in-place, that is, it should modify the `du` array in place.
For details, see the [DifferentialEquations.jl documentation](https://diffeq.sciml.ai/stable/tutorials/sde_example/).

## Example

```julia
using Distributions
c = rand(Community, 3; A_ij=Normal(0, 0.1))
function white_noise!(du, u, p, t)
    for i in eachindex(du)
        du[i] = 0.1 # Noise intensity.
    end
end
u0, tspan = [1.0, 1.0, 1.0], (0, 1_000)
solve(c, u0, tspan, white_noise!)
```
"""
function solve(c::Community, u0, tspan, noise!::Function; kwargs...)
    function f!(du, u, _, _)
        for i in eachindex(u)
            u[i] < 0 && (u[i] = 0)
            du[i] = c.r[i] * u[i] * (1 + sum(c.A[i, :] .* u) / c.K[i])
        end
    end
    prob = SDEProblem(f!, noise!, u0, tspan)
    DifferentialEquations.solve(prob; kwargs...)
end

export solve
