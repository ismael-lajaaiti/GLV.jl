
"""
    cavity_parameters(c::Community)

Take a [`Community`](@ref) and return its summary statistics
for the cavity method. The parameters are:

- `mu`: the mean of the interaction strengths
- `sigma`: the standard deviation of the interaction strengths
- `gamma`: the correlation between the interaction strengths
- `K_std`: the standard deviation of the carrying capacities
- `K_mean`: the mean of the carrying capacities

See also [`cavity_predictions`](@ref).
"""
function cavity_parameters(c::Community)
    A = c.A
    S = richness(c)
    a_ij = [A[i, j] for i in 1:S, j in 1:S if i != j]
    a_ji = [A[j, i] for i in 1:S, j in 1:S if i != j]
    mu = mean(a_ij) * (S - 1)
    sigma = std(a_ij) * sqrt(S - 1)
    gamma = cor(a_ij, a_ji)
    K_mean, K_std = mean(c.K), std(c.K)
    (; mu, sigma, gamma, K_std, K_mean)
end
export cavity_parameters

function cavity_system(x, p)
    mu, sigma, gamma, K_std, K_mean = p
    phi, N_mean, N2_mean, v = x
    u_hat = 1 - phi * gamma * sigma^2 * v
    N0_mean = (K_mean + phi * mu * N_mean) / u_hat
    N0_var = K_std^2 + sigma^2 * phi * N2_mean / u_hat^2
    N0_var <= 1e-5 && (N0_var = 1e-5)
    P(N0) = exp(-(N0 - N0_mean)^2 / N0_var) / sqrt(pi * N0_var)
    [
        phi - quadgk(N0 -> P(N0), 0, Inf)[1],
        N_mean - quadgk(N0 -> N0 * P(N0), 0, Inf)[1] / phi,
        N2_mean - quadgk(N0 -> N0^2 * P(N0), 0, Inf)[1] / phi,
        v - 1 / u_hat
    ]
end


"""
    cavity_predictions(mu, sigma, gamma, K_std; K_mean=1)

Predict the following community properties using cavity method:
- `phi`: the fraction of surviving species
- `N_mean`: the mean of the species abundance distribution
- `N2_mean`: the second moment of the species abundance distribution
- `v`: species response coefficient, which is the derivative of the species abundance with respect to its carrying capacity.

If the solver doesn't converge, all returned values are set to zero.
This usually means that the community is expected to collapse or explode.

# Example

```julia
S = 100
μ, σ = 0, 1
c = rand(Community, S; A_ij=Normal(μ / S, σ / sqrt(S)))
cavity_predictions(c)
```

# References

- [bunin2017](@citet)
- [barbier2017](@citet)
- [barbier2018](@citet)

See also [`cavity_parameters`](@ref).
"""
function cavity_predictions(mu, sigma, gamma, K_std; K_mean=1)
    p = (mu, sigma, gamma, K_std, K_mean)
    u0 = [1, 1, 1, 1]
    prob = NonlinearProblem(cavity_system, u0, p)
    sol = NonlinearSolve.solve(prob)
    if sol.retcode == ReturnCode.Success
        phi, N_mean, N2_mean, v = sol.u
    else
        phi, N_mean, N2_mean, v = 0, 0, 0, 0
    end
    (; phi, N_mean, N2_mean, v)
end

"""
    cavity_predictions(c::Community)

Can take a `Community` and extract the summary statistics for the cavity method.
"""
function cavity_predictions(c::Community)
    mu, sigma, gamma, K_std, K_mean = cavity_parameters(c)
    cavity_predictions(mu, sigma, gamma, K_std; K_mean)
end

export cavity_predictions
