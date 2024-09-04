# # The cavity method
#
# Plot the fraction of surviving species predicted by the cavity methods
# for different values of `mu` and `sigma`.
# ``\\gamma = 1``.

using CairoMakie
using GLV
set_theme!(theme_minimal())

K_std = 0
gamma = 1
n = 10 # Grid size.
sigma_values = LinRange(0.2, 1.2, n) |> collect
mu_values = LinRange(-4, 1.0, n) |> collect
phi_values = [cavity_predictions(mu, sigma, gamma, K_std).phi for mu in mu_values, sigma in sigma_values]

fig = Figure();
ax = Axis(fig[1, 1]; xlabel="μ", ylabel="σ")
heatmap!(mu_values, sigma_values, phi_values)
# ax.yreversed = true
ax.xreversed = true
fig[1, 2] = Colorbar(fig[1, 1])
fig

# When adding `x` and `y` together we obtain a new rational number:

z = x + y
