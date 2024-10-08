"""
    Community(A, r, K)

Create a community with interaction matrix `A`,
growth rates `r`, and carrying capacities `K`.

# Example

```julia
A = [-1 0; 0 -1]
r = [1, 1]
K = [1, 1]
community = Community(A, r, K)
```

See also [`solve`](@ref), [`simulate_pulse`](@ref).
"""
mutable struct Community
    A::AbstractMatrix # Interactions.
    r::AbstractVector # Growth rates.
    K::AbstractVector # Carrying capacities.
    function Community(A, r, K)
        @assert size(A, 1) == size(A, 2) == length(r) == length(K)
        new(A, r, K)
    end
end
export Community

"""
    Base.rand(
    ::Type{Community},
    S::Int;
    A_ij::Distribution=Normal(0, 1),
    r_i::Distribution=Normal(1, 0),
    K_i::Distribution=Normal(1, 0),
    interaction::Symbol=:default,
    )

Genereate a random community with `S` species.
Parameters are drawn from specified distributions.
By default, species growth rates and carrying capacities are set to one.
Species self-regulation, that is the diagonal of `A`, is set to -1.

# Example

Here is a minimal example, where we randomly draw interactions strengths
in a normal distribution.

```julia
using Distributions
c = rand(Community, 10; A_ij = Normal(-1, 0.1))
```

Interactions can also be drawn in a multivariate distribution.
This allows to generate correlated interactions between pair of species.

```julia
using Distributions
S = 100
μ = [-1, -2]
σ = [0.2, 0.1]
ρ = -0.5 # Correlation between A_ij and A_ji.
Σ = [σ[1]^2 ρ*σ[1]*σ[2]; ρ*σ[1]*σ[2] σ[2]^2] # Covariance matrix.
c = rand(Community, S; A_ij = MvNormal(μ, Σ))

# Check that the interaction matrix is correct.
A_ij = [c.A[i, j] for i in 1:S for j in i+1:S]
A_ji = [c.A[j, i] for i in 1:S for j in i+1:S]
mean(A_ij), mean(A_ji)
std(A_ij), std(A_ji)
cor(A_ij, A_ji)
```

See also [`Community`](@ref).
"""
function Base.rand(
    ::Type{Community},
    S::Int;
    A_ij::Distribution = Normal(0, 1),
    r_i::Distribution = Normal(1, 0),
    K_i::Distribution = Normal(1, 0),
    interaction::Symbol = :default,
)
    @assert interaction ∈ [:default, :core]
    r = rand(r_i, S)
    K = rand(K_i, S)
    multivariate_dist = typeof(A_ij) <: MultivariateDistribution
    if !multivariate_dist
        A = rand(A_ij, S, S)
    else
        n = round(Int, S * (S - 1) / 2)
        A_elements = rand(A_ij, n)
        A = zeros(S, S)
        k = 1
        for i in 1:S, j in (i+1):S
            A[i, j] = A_elements[1, k]
            A[j, i] = A_elements[2, k]
            k += 1
        end
    end
    A[diagind(A)] .= -1
    if interaction == :core
        A = Diagonal(K) * A * Diagonal(1 ./ K)
    end
    Community(A, r, K)
end

"""
    abundance(c::Community)

Compute the equilibrium abundance of species in community `c`.
Assumes that `A` is invertible.

## Example

The equilibrium abundance of two non-interacting species
is equal to their carrying capacities.

```jldoctest; output=false
A = [-1 0; 0 -1]
r = [1, 1]
K = [1, 2]
c = Community(A, r, K)
abundance(c) == K

# output

true
```

See also [`relative_yield`](@ref)
"""
abundance(c::Community) = -inv(c.A) * c.K
export abundance

"""
    relative_yield(c::Community)

Compute the equilibrium relative yield of species in community `c`.
Relative yield is the ratio of abundance to carrying capacity.
Assumes that `A` is invertible.

## Example

The equilibrium relative yields of two non-interacting species
are equal to the one.

```jldoctest; output=false
A = [-1 0; 0 -1]
r = [1, 1]
K = [1, 2]
c = Community(A, r, K)
relative_yield(c) == [1, 1]

# output

true
```

See also [`abundance`](@ref).
"""
relative_yield(c::Community) = abundance(c) ./ c.K
export relative_yield

"""
    core_interactions(c::Community)

Compute the 'core' interactions of the community.
Core interactions are the species interactions rescaled
in a relevant manner to study species coexistence.
Formally, the core interactions write

```math
b_{ij} = a_{ij} K_i / K_j
```

where ``a_{ij}`` is the interaction from species ``j`` to species ``i``.

For more information refer to [Barbier and Arnoldi 2017](https://doi.org/10.1101/147728).
"""
function core_interactions(c::Community)
    A, K = c.A, c.K
    Diagonal(1 ./ K) * A * Diagonal(K)
end
export core_interactions

"""
    richness(c::Community)

Species richness of the community `c`.
"""
function richness(c::Community)
    @assert length(c.r) == length(c.K) == size(c.A, 1) == size(c.A, 2)
    length(c.r)
end
export richness

"""
    assemble(c::Community; u0 = ones(richness(c)), tspan = (0, 10_000))

Assemble the pool of species in the community `c`.
Return the subcommunity of species that are alive.
`u0` is the initial condition for the simulation.
`tspan` defines the duration of the simulation.
"""
function assemble(c::Community; u0 = ones(richness(c)), tspan = (0, 10_000))
    sol = solve(c, u0, tspan)
    alive = sol.u[end] .> 1e-6
    Community(c.A[alive, alive], c.r[alive], c.K[alive])
end
export assemble
