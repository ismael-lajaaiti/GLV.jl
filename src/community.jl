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
    r_i::Union{Nothing,Distribution}=nothing,
    K_i::Union{Nothing,Distribution}=nothing,
)

Genereate a random community with `S` species.
Parameters are drawn from specified distributions.
By default, species growth rates and carrying capacities are set to one.
Species self-regulation, that is the diagonal of `A`, is set to -1.

# Example

```julia
using Distributions
c = rand(Community, 10; A_ij=Normal(-1, 0.1))
```

See also [`Community`](@ref).
"""
function Base.rand(
    ::Type{Community},
    S::Int;
    A_ij::Distribution=Normal(0, 1),
    r_i::Distribution=Normal(1, 0),
    K_i::Distribution=Normal(1, 0),
)
    A = rand(A_ij, S, S)
    A[diagind(A)] .= -1
    r = rand(r_i, S)
    K = rand(K_i, S)
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
function abundance(c::Community)
    -inv(c.A) * c.K
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
function relative_yield(c::Community)
    abundance(c) ./ c.K
end
export relative_yield
