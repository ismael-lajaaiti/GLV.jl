```@meta
EditURL = "response-to-pulse.jl"
```

# How biotic interactions structure species responses to pulse perturbations

We will study how interactions between species shape their responses to perturbations,
inspired by this study [lajaaiti2024](@cite).
In particular, we are going to show that the sensitivity of species to perturbations
is strongly related to their relative yield.
Species relative yield is the ratio of the species abundance in the community
over the species abundance when alone (that is, its carrying capacity).
Thus, the relative yield captures the impact of interactions on species equilibrium abundance.

First, let's import the packages we need.

````@example response-to-pulse
using GLV
using LinearAlgebra
using DataFrames
using Distributions
using CairoMakie
set_theme!(theme_minimal())
````

Then, let's create a random competitive community.

````@example response-to-pulse
S = 30
mu, sigma = -1, 0.2
K_std = 0.3
c = rand(
    Community,
    S;
    A_ij = Normal(mu / S, sigma / sqrt(S)),
    K_i = Normal(1, K_std),
    interaction = :core,
)
N_eq = abundance(c)
eta_eq = relative_yield(c)
````

We can check the equilibrium of the community, and ensure that all species
have a positive abundance.

Next, we will perform a series of random pulse perturbations and record each species responses.
But before that, let's define how we quantify species responses to perturbations.

````@example response-to-pulse
"""
    species_responses(sol)

Quantify how strongly each species responds to the simulated perturbation,
whose trajectory is given by `sol`.
The species response is its summed relative deviation from its equilibrium over its recovery.
"""
function species_response(sol, N_eq)
    dist_to_eq = abs.(Array(sol) .- N_eq) ./ N_eq
    mean(dist_to_eq; dims = 2)
end
````

We can now simulate the species responses to pulse perturbations.

````@example response-to-pulse
tspan = (0, 1_000)
n_pulses = 100
df = DataFrame(; species = Int64[], response = Float64[])
for _ in 1:n_pulses
    x = rand(Normal(0, 0.1), S) .* N_eq
    sol = simulate_pulse(c, x, tspan; saveat = 1)
    responses = species_response(sol, N_eq)
    for i in 1:S
        push!(df, (i, responses[i]))
    end
end
````

We then average the species responses over all the pulse realizations,
and add the species abundances and relative yields to the dataframe.
This last step is not necessary, but it will allow us to plot the results.

````@example response-to-pulse
df = combine(groupby(df, :species), :response => mean)
df.abundance = N_eq
df.relative_yield = eta_eq
````

Lastly, we can plot the results.
We see that relative yield, instead of species abundance, strongly control the species responses to pulse perturbations.
Specifically, species with low relative yield are the species that are the most sensitive to pulse perturbations.

````@example response-to-pulse
fig = Figure();
ax1 = Axis(
    fig[1, 1];
    ylabel = "Species response to pulse",
    xlabel = "Species relative yield",
    yscale = log10,
)
scatter!(df.relative_yield, df.response_mean)
ax2 = Axis(fig[1, 2]; xlabel = "Species abundance", yscale = log10)
scatter!(df.abundance, df.response_mean)
fig
````

