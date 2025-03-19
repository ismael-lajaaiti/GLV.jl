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
    u::AbstractVector # Can species grow by themselves.
    θ::AbstractVector # Parameter of θ-logistic.
    function Community(A, r, K, u, θ)
        @assert size(A, 1) == size(A, 2) == length(r) == length(K) == length(u) == length(θ)
        new(A, r, K, u, θ)
    end
end
export Community

"""
    SublinearCommunity(A, r, m, k, B0)

Create a community, governed by the sublinear model.
A is the interaction matrix.
r the vector of species growth rates.
m the vector of species mortality rates.
k the vectory of sublinear scaling.
B0 the vector of lower bound on species densities.
"""
mutable struct SublinearCommunity
    A::AbstractMatrix # Interactions.
    r::AbstractVector # Growth rates.
    m::AbstractVector # Mortality rates.
    k::AbstractVector # Growth scaling.
    B0::AbstractVector # Minimal densities.
    SublinearCommunity(A, r, m, k, B0) = new(A, r, m, k, B0)
end
export SublinearCommunity

function Community(A, r, K)
    u = fill(1, length(r))
    θ = fill(1, length(r))
    Community(A, r, K, u, θ)
end

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
    u = fill(1, S),
    θ = fill(1, S),
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
    Community(A, r, K, u, θ)
end

"""
    Base.rand(
    ::Type{SublinearCommunity},
    S::Int;
    A_ij::Distribution = Normal(0, 1),
    r_i::Distribution = Normal(1, 0),
    m_i::Distribution = Normal(0.1, 0),
    k_i::Distribution = Normal(0.75, 0),
    B0_i::Distribution = Normal(0.01, 0),
    interaction::Symbol = :default,

)

Generate a [`SublinearCommunity`](@ref), in which species follow a sublinear growth,
with random parameters.
"""
function Base.rand(
    ::Type{SublinearCommunity},
    S::Int;
    A_ij::Distribution = Normal(0, 1),
    r_i::Distribution = Normal(1, 0),
    m_i::Distribution = Normal(10^(-0.5), 0),
    k_i::Distribution = Normal(0.75, 0),
    B0_i::Distribution = Normal(0.01, 0),
    interaction::Symbol = :default,
)
    @assert interaction ∈ [:default, :core]
    r = rand(r_i, S)
    m = rand(m_i, S)
    k = rand(k_i, S)
    B0 = rand(B0_i, S)
    K = carrying_capacity(r, m, k, B0)
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
    A[diagind(A)] .= 0
    if interaction == :core
        A = Diagonal(K) * A * Diagonal(1 ./ K)
    end
    SublinearCommunity(A, r, m, k, B0)
end

"""
    carrying_capacity(c::SublinearCommunity)

Compute species carrying capacities, that is, their abundance when alone.
"""
carrying_capacity(c::SublinearCommunity) = carrying_capacity(c.r, c.m, c.k, c.B0)
carrying_capacity(r, m, k, B0) = B0 .* (r ./ m) .^ (1 ./ (1 .- k))
export carrying_capacity

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
function abundance(c::Community; t_end = 100_000)
    logistic = all(c.θ .== 1)
    if logistic # Can be solved analytically.
        Beq = -inv(c.A) * (c.u .* c.K)
    else
        Beq = solve(c, c.K, (0, t_end))[end]
    end
    Beq
end
function abundance(c::SublinearCommunity)
    K = carrying_capacity(c)
    solve(c, K, (0, 100_000))[end]
end
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
relative_yield(c::Community) = abundance(c) ./ (c.K .* c.u)
relative_yield(c::SublinearCommunity) = abundance(c) ./ carrying_capacity(c)
export relative_yield

relative_selfregulation(c::Community) = (abundance(c) ./ c.K) .^ c.θ
function relative_selfregulation(c::SublinearCommunity)
    B = abundance(c)
    B_dfdB = (1 .- c.k) .* c.r .* c.B0 .^ (1 .- c.k) .* B .^ (c.k .- 1)
    f0 = c.r .- c.m
    B_dfdB ./ f0
end
export relative_selfregulation

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
function core_interactions(c::SublinearCommunity)
    A = c.A
    K = carrying_capacity(c)
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
