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
    f!(du, u, c, _) =
        for i in eachindex(u)
            u[i] < 0 && (u[i] = 0)
            if u[i] < 1e-12
                du[i] = 0
            else
                du[i] =
                    c.r[i] *
                    u[i] *
                    (c.u[i] + (u[i]^(c.θ[i] - 1) / c.K[i]^c.θ[i]) * sum(c.A[i, :] .* u))
            end
        end
    prob = ODEProblem(f!, u0, tspan, c)
    DifferentialEquations.solve(prob; kwargs...)
end

function solve(c::SublinearCommunity, u0, tspan; kwargs...)
    f!(du, u, c, _) =
        for i in eachindex(u)
            u[i] < 0 && (u[i] = 0) # Species cannot have negative abundances.
            du[i] =
                u[i] * (
                    c.r[i] * (c.B0[i] / u[i])^(1 - c.k[i]) + # Sublinear growth.
                    -c.m[i] + # Mortality.
                    sum(c.A[i, :] .* u) # Interactions.
                )
            u[i] < c.B0[i] && (du[i] = 0) # Lower bound on growth.
        end
    prob = ODEProblem(f!, u0, tspan, c)
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
c = rand(Community, 3; A_ij = Normal(0, 0.1))
white_noise!(du, u, p, t) =
    for i in eachindex(du)
        du[i] = 0.1 # Noise intensity.
    end
u0, tspan = [1.0, 1.0, 1.0], (0, 1_000)
solve(c, u0, tspan, white_noise!)
```
"""
function solve(c::Community, u0, tspan, noise!::Function; kwargs...)
    f!(du, u, c, _) =
        for i in eachindex(u)
            u[i] < 0 && (u[i] = 0)
            du[i] =
                c.r[i] *
                u[i] *
                (c.u[i] + (u[i]^(c.θ[i] - 1) / c.K[i]^c.θ[i]) * sum(c.A[i, :] .* u))
        end
    prob = SDEProblem(f!, noise!, u0, tspan, c)
    DifferentialEquations.solve(prob; kwargs...)
end

export solve
