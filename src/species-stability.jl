"""
    species_reactivity(c::Community)

Compute the species reactivity for each species in the community `c`.
Species reactivity correspond to the worst initial response (given by the slope)
to a pulse perturbation.
Formally, it writes

```math
R_0^{(i)} = \\sqrt{\\sum_{j\\neq i} a_{ij}^2 \\eta_j^2}
```

where ``a_{ij}`` is the interaction from species ``j`` to species ``i``
and ``\\eta`` is the relative yield.

## Example

```julia
using Distributions
S = 50
μ, σ = -1, 0.2
c = rand(
    Community,
    S;
    A_ij = Normal(μ / S, σ / sqrt(S)),
    K_i = Uniform(1, 10),
    interaction = :core,
)
r0 = species_reactivity(c)
cor(abundance(c), r0) # Weak negative correlation.
cor(relative_yield(c), r0) # Strong negative correlation.
```

For more information refer to [Lajaaiti et al. 2024](https://doi.org/10.22541/au.172479485.55602895/v1).
#TODO: Update reference when the article is in press.
"""
function species_reactivity(c::Community)
    S = size(c.A, 1)
    η = relative_yield(c)
    B = core_interactions(c)
    B_nodiag = B - Diagonal(B)
    [sqrt(sum(B_nodiag[i, :] .* η) .^ 2) for i in 1:S]
end
export species_reactivity

"""
    get_species_cv(sol::AbstractODESolution)

Compute the coefficient of variation of each species given stochastic trajectories.

# Example

```julia
S = 10
c = rand(Community, S; A_ij = Normal(0, 0))
N_eq = abundance(c)
noise! = function (du, u, p, t) for i in eachindex(u)
    du[i] = u[i] * 0.1
    end
end
tspan = (0, 100)
sol = solve(c, N_eq, tspan, noise!)
get_species_cv(sol)
```

See also [`solve`](@ref).
"""
get_species_cv(sol) = vec(std(sol; dims = 2) ./ mean(sol; dims = 2))
export get_species_cv
