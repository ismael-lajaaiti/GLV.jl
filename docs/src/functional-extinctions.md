```@meta
EditURL = "functional-extinctions.jl"
```

# Functional extinctions in ecological communities

We can distinguish two types of species extinctions.
First, numerical extinctions occur when the very last member of a species dies.
Second, functional extinctions occur when the species is still present
in the community, but is too rare to fulfill its function.

Following methods of [saterberg2013](@cite) we consider that a species goes
functionall extinct when a decrease in its abundance results in the extinction
of some other species in the community.
Here a will consider plant communities, thus to decrease species abundances
we decrease their carrying capacity.

````@example functional-extinctions
using GLV
using LinearAlgebra
using Distributions
using CairoMakie
using Random
set_theme!(theme_minimal())
Random.seed!(123)


S = 30
mu, sigma = -1 / S, 0.2 / sqrt(S)
K_mean, K_std = 1, 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu, sigma),
    K_i = Normal(K_mean, K_std),
    interaction = :core,
)
````

We can check that the community results in a stable equilibrium,
by checking that all species have positive abundances.

````@example functional-extinctions
all(abundance(c) .> 0)
````

Now let's compute the ecologically effective population size (EEP) of each species,
that is, the population size below which some other species in the community go extinct.
For each species $j$, we want to find the smallest decrease in its growth rate
that would result in the extinction of some other species $i$.
It can be shown that this quantity writes

```math
\varepsilon_j = \min_i(-\frac{\hat{N}_i}{A^{-1}_{ij}} | \varepsilon_j > 0)
```

Let's compute this quantity for species each species.

````@example functional-extinctions
Neq = abundance(c)
epsilon = zeros(S)
for j in 1:S
    eps = -abundance(c) ./ inv(c.A)[:, j]
    epsilon[j] = minimum(eps[eps.>0])
end
epsilon
````

We now have the smallest derease in carrying capacity that would result in the extinction
of a species. The species becoming extinct can be the focal species or some other species.
In the first scenario, there is no functional extinction but solely a numerical extinction.

We can compute the ecologically effective population size (EEP) of each species, that is,
the minimal population size below which some species go extinct.
When there is no functional extinction, but a numerical extinction,
the EEP size is equal to zero.

````@example functional-extinctions
delta_N = Diagonal(inv(c.A)) * epsilon
N_EEP = Neq .+ delta_N
````

We see that all species EEP size are zero.
This is not suprising because we have set very weak interactions strengths.
Thus, the species are not dependent on each other, and the decrease in abundance
of one species has very little effect on the others.

Let's repeat our analysis with a stronger interaction matrix.

````@example functional-extinctions
mu, sigma = -1 / S, 0.8 / sqrt(S)
K_mean, K_std = 1, 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu, sigma),
    K_i = Normal(K_mean, K_std),
    interaction = :core,
)
````

Because interactions are stronger, we expect some extinctions during the community assembly.

````@example functional-extinctions
Neq = abundance(c)
all(Neq .> 0)
````

So let's assemble the community, and keep only the surviving species.
This can be done simply by calling the [`assemble`](@ref) function
that is been designed precisely for this purpose.

````@example functional-extinctions
c_new = assemble(c)
````

We can check that some species have gone extinct during the assembly.

````@example functional-extinctions
S_new = richness(c_new) # Smaller than the initial richness S.
````

````@example functional-extinctions
epsilon = zeros(S_new)
for j in 1:S_new
    eps = -abundance(c_new) ./ inv(c_new.A)[:, j]
    epsilon[j] = minimum(eps[eps.>0])
end
delta_N = Diagonal(inv(c_new.A)) * epsilon
Neq = abundance(c_new)
N_EEP = Neq .+ delta_N
````

We see that all EEP sizes are not zeros.
This means that some species are functionally necessary for the communtiy.

Finally, we can plot the EEP size of each species against their initial equilibrium abundance.

````@example functional-extinctions
fig = Figure();
ax = Axis(fig[1, 1]; xlabel = "Equilibrium abundance", ylabel = "EEP size")
positive_eep = N_EEP .> 0
scatter!(Neq[positive_eep], N_EEP[positive_eep])
fig
````

We note a stong positive trend, that is, the more abundant the species the larger its EEP size
therefore the more "quickly" it goes functionally extinct.

