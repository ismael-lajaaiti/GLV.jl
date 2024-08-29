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

See also [`solve`](@ref).
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
