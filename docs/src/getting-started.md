```@meta
CurrentModule = GLV
```

# Getting started

## Installation

The package can be installed with

```julia
using Pkg
Pkg.add("GLV")
```

Then the library can be loaded with

```@example main
using GLV
```


## Create a community

The first step to simulate the dynamics of a GLV community
is to define its parameters.
The community parameters are the following:
- `A` the interaction matrix
- `r` the vector of species growth rates
- `K` the vector of species carrying capacities

Therefore a simple community can be created with

```@example main
A = [-1 0.1; -0.1 -1]
r = [0.1, 0.2]
K = [2, 3]
c = Community(A, r, K)
```

Yet, because GLV communities are often large and their
parameters drawn from random distributions, we provide
a [`rand`](@ref) method to create random communities.


```@example main
using Distributions # To draw random numbers.
S = 10 # Number of species.
A_ij = Normal(0, 0.1) # Distribution of interaction coefficients.
c = rand(Community, S; A_ij = A_ij)
c.A
```

We can see that interspecific interactions have been drawn in the given distribution,
and that intraspecific interactions are set to -1.

Moreover, by default, species growth rates and carrying capacities are set 1.

```@example main
c.r, c.K
```



## Simulate the community dynamics

