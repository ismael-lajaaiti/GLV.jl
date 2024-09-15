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

The Lotka-Volterra model writes

```math
\frac{dN_i}{dt} = r_i N_i \left(1 + \frac{\sum_j A_{ij} N_j} {K_i}\right)
```

where $N_i$ is the abundance of species $i$, $r_i$ its growth rate,
and $K_i$ its carrying capacity.

Therefore a [`Community`](@ref) is an object that stores the following
parameters:
- `A` the interaction matrix
- `r` the vector of species growth rates
- `K` the vector of species carrying capacities

For example, to create a simple community

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

## Simulate its dynamics

Once the community is created, its dynamics can be simulated using
the `solve` function

```@example main
N0 = rand(S) # Initial abundances.
tspan = (0, 50) # Simulation duration.
sol = solve(c, N0, tspan)
```

These trajectories can be plotted using your favourite plotting library.
For example, with CairoMakie

```@example main
using CairoMakie
lines(sol)
```

The species equilibrium abundances can also be computed analytically with

```@example main
abundance(c)
```

!!! note
    Equilibrium abundances are computed by solving the system of ODEs
    when the time derivative of the abundances is zero.
    That is, $N^* = - A^{-1} \mathbf{K}$.

We can check that these values correspond to the abundances
obtained at the end of the simulation

```@example main
sol[end]
```

For more information and advanced usage see the examples.
