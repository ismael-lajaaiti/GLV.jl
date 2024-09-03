
"""
    cavity_parameters(c::Community)

TBW
"""
function cavity_parameters(c::Community)
    A = c.A
    S = richness(c)
    a_ij = [A[i, j] for i in 1:S, j in 1:S if i != j]
    a_ji = [A[j, i] for i in 1:S, j in 1:S if i != j]
    mu = mean(a_ij) * (S - 1)
    sigma = std(a_ij) * sqrt(S - 1)
    gamma = cor(a_ij, a_ji)
    (; mu, sigma, gamma)
end
export cavity_parameters

function cavity_system(x, c::Community)
    # Parameters that we want to solve for.
    phi = x[1]
    N_mean = x[2]
    N2_mean = x[3]
    v = x[4]
    # Parameters defined for convenience.
    S = richness(c)
    mu, sigma, gamma = cavity_parameters(c)
    u_hat = 1 + mu / S - phi * gamma * sigma^2 * v
    K_mean, K_std = mean(c.K), std(c.K)
    N0_mean = (K_mean + phi * mu * N_mean) / u_hat
    N0_var = K_std^2 + sigma^2 * phi * N2_mean / u_hat^2
    N0_var <= 1e-5 && (N0_var = 1e-5)
    P(N0) = exp(-(N0 - N0_mean)^2 / N0_var) / sqrt(pi * N0_var)
    # @info N0_mean, N0_var, u_hat
    # System to solve.
    [
        phi - quadgk(N0 -> P(N0), 0, Inf)[1],
        N_mean - quadgk(N0 -> N0 * P(N0), 0, Inf)[1] / phi,
        N2_mean - quadgk(N0 -> N0^2 * P(N0), 0, Inf)[1] / phi,
        v - 1 / u_hat
    ]
end

"""
    cavity_predictions(c::Community)

TBW

## Example

```julia

using Distributions

S = 100
μ, σ = -3, 1
c = rand(Community, S; A_ij=Normal(μ / S, σ / sqrt(S)))
cavity_predictions(c)

N0, tspan = fill(1, S), (0, 1_000)
u = GLV.solve(c, N0, tspan)
u_end = u[end]
mean(u[end] .> 0)
mean(u_end[u_end .> 0])

```
"""
function cavity_predictions(c::Community)
    u0 = [1, 1, 1, 1]
    prob = NonlinearProblem(cavity_system, u0, c)
    NonlinearSolve.solve(prob).u
end
export cavity_predictions
