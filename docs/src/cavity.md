```@meta
EditURL = "cavity.jl"
```

# The cavity method

The cavity method is a mathematical technique originally developed in
statistical physics to study disordered systems like spin glasses.
In recent years, it has been adapted to ecological models to
analyze complex ecosystems with many interacting species.

To illustrate the strength of the cavity method, we will use it to
predict the first two moments of the species abundance
distribution in the assembled community.

First, we have to import few packages.

````@example cavity
using GLV
using DataFrames
using Distributions
using CairoMakie
set_theme!(theme_minimal())
````

`DataFrames` is used to store generated data.
`Distributions` is used to draw community parameters.
`CairoMakie` is our library of choice for plotting.

We set the model parameters.
We will vary the mean interaction strength ``\mu``.

````@example cavity
K_std = 0.2
gamma = 0 # Uncorrelated interactions.
sigma = 0.3
mu_values = -collect(LinRange(0, 2, 50))
````

We then define the parameters of the simulation.

````@example cavity
S = 50
N0 = fill(1, S)
tspan = (0, 1000)
````

Then we can run our simulation, and compare the results with the cavity method.

````@example cavity
df = DataFrame(; simulation=Float64[], prediction=Float64[], parameter=Symbol[])
for mu in mu_values
    c = rand(Community, S; A_ij=Normal(mu / S, sigma / sqrt(S)), K_i=Normal(1, K_std))
    sol = solve(c, N0, tspan)
    N = sol.u[end][sol.u[end].>0]
    N_mean = mean(N)
    N2_mean = mean(N .^ 2)
    simulation = (; N_mean, N2_mean)
    p = cavity_predictions(c)
    prediction = (N_mean=p.N_mean * p.phi, N2_mean=p.N2_mean * p.phi)
    for param in [:N_mean, :N2_mean]
        push!(df, (simulation[param], prediction[param], param))
    end
end
````

And here are the results.

````@example cavity
fig = Figure(size=(500, 300));
title_list = (N_mean="Mean abundance", N2_mean="Mean abundance squared")
for (i, param) in enumerate(unique(df.parameter))
    ylabel = i == 1 ? "Prediction" : ""
    title = title_list[param]
    ax = Axis(fig[1, i]; xlabel="Simulation", ylabel, title)
    scatter!(df[df.parameter.==param, :simulation], df[df.parameter.==param, :prediction])
    ablines!(0, 1; color=:black)
end
fig
````

