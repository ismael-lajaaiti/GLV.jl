var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Barbier, M. and Arnoldi, J.-F. (2017). The cavity method for community ecology, bioRxiv, 147728.\n\n\n\nBarbier, M.; Arnoldi, J.-F.; Bunin, G. and Loreau, M. (2018). Generic assembly patterns in complex ecological communities. Proceedings of the National Academy of Sciences 115, 2156–2161.\n\n\n\nBunin, G. (2017). Ecological communities with Lotka-Volterra dynamics. Phys. Rev. E 95, 042414.\n\n\n\nLajaaiti, I.; Kefi, S. and Arnoldi, J.-F. (2024). How biotic interactions structure species' responses to perturbations.\n\n\n\nSäterberg, T.; Sellman, S. and Ebenman, B. (2013). High frequency of functional extinctions in ecological networks. Nature 499, 468–470.\n\n\n\n","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"EditURL = \"cavity.jl\"","category":"page"},{"location":"cavity/#The-cavity-method","page":"The cavity method","title":"The cavity method","text":"","category":"section"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"The cavity method is a mathematical technique originally developed in statistical physics to study disordered systems like spin glasses. In recent years, it has been adapted to ecological models to analyze complex ecosystems with many interacting species.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"To illustrate the strength of the cavity method, we will use it to predict the first two moments of the species abundance distribution in the assembled community.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"First, we have to import few packages.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"using GLV\nusing DataFrames\nusing Distributions\nusing CairoMakie\nset_theme!(theme_minimal())","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"DataFrames is used to store generated data. Distributions is used to draw community parameters. CairoMakie is our library of choice for plotting.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"We set the model parameters. We will vary the mean interaction strength mu.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"K_std = 0.2\ngamma = 0 # Uncorrelated interactions.\nsigma = 0.3\nmu_values = -collect(LinRange(0, 2, 50))","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"We then define the parameters of the simulation.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"S = 50\nN0 = fill(1, S)\ntspan = (0, 1000)","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"Then we can run our simulations and store the generated data.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"df = DataFrame(; simulation=Float64[], prediction=Float64[], parameter=Symbol[])\nfor mu in mu_values\n    c = rand(Community, S; A_ij=Normal(mu / S, sigma / sqrt(S)), K_i=Normal(1, K_std))\n    sol = solve(c, N0, tspan)\n    N = sol.u[end][sol.u[end].>0]\n    N_mean = mean(N)\n    N2_mean = mean(N .^ 2)\n    simulation = (; N_mean, N2_mean)\n    p = cavity_predictions(c)\n    prediction = (N_mean=p.N_mean * p.phi, N2_mean=p.N2_mean * p.phi)\n    for param in [:N_mean, :N2_mean]\n        push!(df, (simulation[param], prediction[param], param))\n    end\nend","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"Finally, we can plot the results.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"fig = Figure(size=(500, 300));\ntitle_list = (N_mean=\"Mean abundance\", N2_mean=\"Mean abundance squared\")\nfor (i, param) in enumerate(unique(df.parameter))\n    ylabel = i == 1 ? \"Prediction\" : \"\"\n    title = title_list[param]\n    ax = Axis(fig[1, i]; xlabel=\"Simulation\", ylabel, title)\n    scatter!(df[df.parameter.==param, :simulation], df[df.parameter.==param, :prediction])\n    ablines!(0, 1; color=:black)\nend\nfig","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"We can that the predictions are very close to the simulations.","category":"page"},{"location":"cavity/","page":"The cavity method","title":"The cavity method","text":"For more details about the cavity method, you can check these References.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"EditURL = \"response-to-pulse.jl\"","category":"page"},{"location":"response-to-pulse/#How-biotic-interactions-structure-species-responses-to-pulse-perturbations","page":"Interactions structure species responses","title":"How biotic interactions structure species responses to pulse perturbations","text":"","category":"section"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"We will study how interactions between species shape their responses to perturbations, inspired by this study (Lajaaiti et al., 2024). In particular, we are going to show that the sensitivity of species to perturbations is strongly related to their relative yield. Species relative yield is the ratio of the species abundance in the community over the species abundance when alone (that is, its carrying capacity). Thus, the relative yield captures the impact of interactions on species equilibrium abundance.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"First, let's import the packages we need.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"using GLV\nusing LinearAlgebra\nusing DataFrames\nusing Distributions\nusing CairoMakie\nset_theme!(theme_minimal())","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"Then, let's create a random competitive community.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"S = 30\nmu, sigma = -1, 0.2\nK_std = 0.3\nc = rand(\n    Community,\n    S;\n    A_ij = Normal(mu / S, sigma / sqrt(S)),\n    K_i = Normal(1, K_std),\n    interaction = :core,\n)\nN_eq = abundance(c)\neta_eq = relative_yield(c)","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"We can check the equilibrium of the community, and ensure that all species have a positive abundance.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"Next, we will perform a series of random pulse perturbations and record each species responses. But before that, let's define how we quantify species responses to perturbations.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"\"\"\"\n    species_responses(sol)\n\nQuantify how strongly each species responds to the simulated perturbation,\nwhose trajectory is given by `sol`.\nThe species response is its summed relative deviation from its equilibrium over its recovery.\n\"\"\"\nfunction species_response(sol, N_eq)\n    dist_to_eq = abs.(Array(sol) .- N_eq) ./ N_eq\n    mean(dist_to_eq; dims = 2)\nend","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"We can now simulate the species responses to pulse perturbations.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"tspan = (0, 1_000)\nn_pulses = 100\ndf = DataFrame(; species = Int64[], response = Float64[])\nfor _ in 1:n_pulses\n    x = rand(Normal(0, 0.1), S) .* N_eq\n    sol = simulate_pulse(c, x, tspan; saveat = 1)\n    responses = species_response(sol, N_eq)\n    for i in 1:S\n        push!(df, (i, responses[i]))\n    end\nend","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"We then average the species responses over all the pulse realizations, and add the species abundances and relative yields to the dataframe. This last step is not necessary, but it will allow us to plot the results.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"df = combine(groupby(df, :species), :response => mean)\ndf.abundance = N_eq\ndf.relative_yield = eta_eq","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"Lastly, we can plot the results. We see that relative yield, instead of species abundance, strongly control the species responses to pulse perturbations. Specifically, species with low relative yield are the species that are the most sensitive to pulse perturbations.","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"fig = Figure();\nax1 = Axis(\n    fig[1, 1];\n    ylabel = \"Species response to pulse\",\n    xlabel = \"Species relative yield\",\n    yscale = log10,\n)\nscatter!(df.relative_yield, df.response_mean)\nax2 = Axis(fig[1, 2]; xlabel = \"Species abundance\", yscale = log10)\nscatter!(df.abundance, df.response_mean)\nfig","category":"page"},{"location":"response-to-pulse/","page":"Interactions structure species responses","title":"Interactions structure species responses","text":"If you want to learn more about this, I invite you to read the study (Lajaaiti et al., 2024).","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GLV","category":"page"},{"location":"#Generalized-Lotka-Volterra-Model-in-Julia","page":"Home","title":"Generalized Lotka-Volterra Model in Julia","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for GLV.jl.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed with","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg; Pkg.add(\"GLV\")","category":"page"},{"location":"#Docstrings","page":"Home","title":"Docstrings","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [GLV]","category":"page"},{"location":"#GLV.Community","page":"Home","title":"GLV.Community","text":"Community(A, r, K)\n\nCreate a community with interaction matrix A, growth rates r, and carrying capacities K.\n\nExample\n\nA = [-1 0; 0 -1]\nr = [1, 1]\nK = [1, 1]\ncommunity = Community(A, r, K)\n\nSee also solve, simulate_pulse.\n\n\n\n\n\n","category":"type"},{"location":"#Base.rand-Tuple{Type{Community}, Int64}","page":"Home","title":"Base.rand","text":"Base.rand(\n::Type{Community},\nS::Int;\nA_ij::Distribution=Normal(0, 1),\nr_i::Distribution=Normal(1, 0),\nK_i::Distribution=Normal(1, 0),\ninteraction::Symbol=:default,\n)\n\nGenereate a random community with S species. Parameters are drawn from specified distributions. By default, species growth rates and carrying capacities are set to one. Species self-regulation, that is the diagonal of A, is set to -1.\n\nExample\n\nusing Distributions\nc = rand(Community, 10; A_ij = Normal(-1, 0.1))\n\nSee also Community.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.abundance-Tuple{Community}","page":"Home","title":"GLV.abundance","text":"abundance(c::Community)\n\nCompute the equilibrium abundance of species in community c. Assumes that A is invertible.\n\nExample\n\nThe equilibrium abundance of two non-interacting species is equal to their carrying capacities.\n\nA = [-1 0; 0 -1]\nr = [1, 1]\nK = [1, 2]\nc = Community(A, r, K)\nabundance(c) == K\n\nSee also relative_yield\n\n\n\n\n\n","category":"method"},{"location":"#GLV.assemble-Tuple{Community}","page":"Home","title":"GLV.assemble","text":"assemble(c::Community; u0 = ones(richness(c)), tspan = (0, 10_000))\n\nAssemble the pool of species in the community c. Return the subcommunity of species that are alive. u0 is the initial condition for the simulation. tspan defines the duration of the simulation.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.cavity_parameters-Tuple{Community}","page":"Home","title":"GLV.cavity_parameters","text":"cavity_parameters(c::Community)\n\nTake a Community and return its summary statistics for the cavity method. The parameters are:\n\nmu: the mean of the interaction strengths\nsigma: the standard deviation of the interaction strengths\ngamma: the correlation between the interaction strengths\nK_std: the standard deviation of the carrying capacities\nK_mean: the mean of the carrying capacities\n\nSee also cavity_predictions.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.cavity_predictions-NTuple{4, Any}","page":"Home","title":"GLV.cavity_predictions","text":"cavity_predictions(mu, sigma, gamma, K_std; K_mean=1, kwargs...)\n\nPredict the following community properties using cavity method:\n\nphi: the fraction of surviving species\nN_mean: the mean of the species abundance distribution\nN2_mean: the second moment of the species abundance distribution\nv: species response coefficient, which is the derivative of the species abundance with respect to its carrying capacity.\n\nIf the solver doesn't converge, all returned values are set to zero. This usually means that the community is expected to collapse or explode.\n\nKeyword arguments of NonlinearSolve.solve can be directly passed to this function.\n\nExample\n\nS = 100\nμ, σ = 0, 1\nc = rand(Community, S; A_ij=Normal(μ / S, σ / sqrt(S)))\ncavity_predictions(c)\n\nReferences\n\nBunin (2017)\nBarbier and Arnoldi (2017)\nBarbier et al. (2018)\n\nSee also cavity_parameters.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.cavity_predictions-Tuple{Community}","page":"Home","title":"GLV.cavity_predictions","text":"cavity_predictions(c::Community)\n\nCan take a Community and extract the summary statistics for the cavity method.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.core_interactions-Tuple{Community}","page":"Home","title":"GLV.core_interactions","text":"core_interactions(c::Community)\n\nCompute the 'core' interactions of the community. Core interactions are the species interactions rescaled in a relevant manner to study species coexistence. Formally, the core interactions write\n\nb_ij = a_ij K_i  K_j\n\nwhere a_ij is the interaction from species j to species i.\n\nFor more information refer to Barbier and Arnoldi 2017.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.offdiag-Tuple{Any}","page":"Home","title":"GLV.offdiag","text":"offdiag(A)\n\nDictionnary of off-diagonal elements of a matrix A. Keys are pairs of indices (i, j) and values are the corresponding elements of A.\n\nA = [1 2; 3 4]\noffdiag(A)\n\n\n\n\n\n","category":"method"},{"location":"#GLV.relative_yield-Tuple{Community}","page":"Home","title":"GLV.relative_yield","text":"relative_yield(c::Community)\n\nCompute the equilibrium relative yield of species in community c. Relative yield is the ratio of abundance to carrying capacity. Assumes that A is invertible.\n\nExample\n\nThe equilibrium relative yields of two non-interacting species are equal to the one.\n\nA = [-1 0; 0 -1]\nr = [1, 1]\nK = [1, 2]\nc = Community(A, r, K)\nrelative_yield(c) == [1, 1]\n\nSee also abundance.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.richness-Tuple{Community}","page":"Home","title":"GLV.richness","text":"richness(c::Community)\n\nSpecies richness of the community c.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.simulate_extinctions-Tuple{Community, Any, Any}","page":"Home","title":"GLV.simulate_extinctions","text":"simulate_extinctions(c::Community, idx, tspan)\n\nSimulate the dynamics of the community c after the extinction of the species of indices idx for the time span tspan.\n\nExample\n\nusing Distributions\nc = rand(Community, 5; A_ij = Normal(0, 0.1))\nsimulate_extinctions(c, [1, 3], (0, 100)) # Species 1 and 3 go extinct.\n\nSee also solve, simulate_pulse.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.simulate_noise-Tuple{Community, Function, Any}","page":"Home","title":"GLV.simulate_noise","text":"simulate_noise(c::Community, noise!::Function, tspan)\n\nSimulate the dynamics of the community c with stochastic noise, around its equilibrium. The species equilibrium abundances are given by abundance.\n\nSee also solve, simulate_pulse, simulate_extinctions, simulate_press.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.simulate_press-Tuple{Community, Any, Any}","page":"Home","title":"GLV.simulate_press","text":"simulate_press(c::Community, K_new, tspan; kwargs...)\n\nSimulate the dynamics of the community c following a press perturbation The press perturbation is modeled by a change in the carrying capacities. The species carrying capacities after the perturbation are given by K_new.\n\nExample\n\nusing Distributions\nc = rand(Community, 5; A_ij = Normal(0, 0.1))\nK_new = c.K .- [0.9, 0, 0, 0, 0] # Lower the carrying capacity of the first species.\nsimulate_press(c, K_new, (0, 100))\n\nSee also solve, simulate_pulse, simulate_extinctions.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.simulate_pulse-Tuple{Community, Any, Any}","page":"Home","title":"GLV.simulate_pulse","text":"simulate_pulse(c::Community, x, tspan)\n\nSimulate the recovery of the community c after the pulse perturbation x for the time span tspan. In other words, it simulate the community dynamics of initial conditions N^* + x where N^* is the vector of species equilibrium abundances.\n\nExample\n\nusing Distributions\nS = 5 # Number of species.\nc = rand(Community, S; A_ij = Normal(0, 0.1))\nx = rand(Normal(-2, 0.1), S)\nu = simulate_pulse(c, x, (0, 100))\n\nSee also solve.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.solve-Tuple{Community, Any, Any, Function}","page":"Home","title":"GLV.solve","text":"solve(c::Community, u0, tspan, noise!::Function; kwargs...)\n\nRun solve with stochastic noise, given by the function noise. The noise! function should be defined as in-place, that is, it should modify the du array in place. For details, see the DifferentialEquations.jl documentation.\n\nExample\n\nusing Distributions\nc = rand(Community, 3; A_ij=Normal(0, 0.1))\nfunction white_noise!(du, u, p, t)\n    for i in eachindex(du)\n        du[i] = 0.1 # Noise intensity.\n    end\nend\nu0, tspan = [1.0, 1.0, 1.0], (0, 1_000)\nsolve(c, u0, tspan, white_noise!)\n\n\n\n\n\n","category":"method"},{"location":"#GLV.solve-Tuple{Community, Any, Any}","page":"Home","title":"GLV.solve","text":"DifferentialEquations.solve(c::Community, u0, tspan; kwargs...)\n\nRun the GLV model for community c with initial conditions u0 and time span tspan. The GLV model writes\n\nfracmathrmd N_imathrmdt = r_i N_i left(fracsum_jneq i A_ij N_j - N_iK_iright)\n\nwhere r is the growth rate, A is the interaction matrix, and K is the carrying capacity.\n\nExample\n\nTwo non-interacting species with different carrying capacities.\n\nA = [-1 0; 0 -1] # Only self-interactions.\nr = [1.0, 1.0]\nK = [1.0, 2.0]\nc = Community(A, r, K)\nu0, tspan = [1.0, 1.0], (0, 10_000) # Simulation parameters.\nsol = solve(c, u0, tspan) # Simulate the dynamics.\n\nSee also Community.\n\n\n\n\n\n","category":"method"},{"location":"#GLV.species_reactivity-Tuple{Community}","page":"Home","title":"GLV.species_reactivity","text":"species_reactivity(c::Community)\n\nCompute the species reactivity for each species in the community c. Species reactivity correspond to the worst initial response (given by the slope) to a pulse perturbation. Formally, it writes\n\nR_0^(i) = sqrtsum_jneq i a_ij^2 eta_j^2\n\nwhere a_ij is the interaction from species j to species i and eta is the relative yield.\n\nExample\n\nusing Distributions\nS = 50\nμ, σ = -1, 0.2\nc = rand(Community, S; A_ij=Normal(μ / S, σ / sqrt(S)), K_i=Uniform(1, 10), interaction=:core)\nr0 = species_reactivity(c)\ncor(abundance(c), r0) # Weak negative correlation.\ncor(relative_yield(c), r0) # Strong negative correlation.\n\nFor more information refer to Lajaaiti et al. 2024. #TODO: Update reference when the article is in press.\n\n\n\n\n\n","category":"method"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"EditURL = \"functional-extinctions.jl\"","category":"page"},{"location":"functional-extinctions/#Functional-extinctions-in-ecological-communities","page":"Functional extinctions","title":"Functional extinctions in ecological communities","text":"","category":"section"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We can distinguish two types of species extinctions. First, numerical extinctions occur when the very last member of a species dies. Second, functional extinctions occur when the species is still present in the community, but is too rare to fulfill its function.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Following methods of (Säterberg et al., 2013) we consider that a species goes functionall extinct when a decrease in its abundance results in the extinction of some other species in the community. Here a will consider plant communities, thus to decrease species abundances we decrease their carrying capacity.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"using GLV\nusing LinearAlgebra\nusing Distributions\nusing CairoMakie\nset_theme!(theme_minimal())\n\nS = 30\nmu, sigma = -1 / S, 0.2 / sqrt(S)\nK_mean, K_std = 1, 0.3\nc = rand(\n    Community,\n    S;\n    A_ij = Normal(mu, sigma),\n    K_i = Normal(K_mean, K_std),\n    interaction = :core,\n)","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We can check that the community results in a stable equilibrium, by checking that all species have positive abundances.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"all(abundance(c) .> 0)","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Now let's compute the ecologically effective population size (EEP) of each species, that is, the population size below which some other species in the community go extinct. For each species j, we want to find the smallest decrease in its growth rate that would result in the extinction of some other species i. It can be shown that this quantity writes","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"varepsilon_j = min_i(-frachatN_iA^-1_ij  varepsilon_j  0)","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Let's compute this quantity for species each species.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Neq = abundance(c)\nepsilon = zeros(S)\nfor j = 1:S\n    eps = -abundance(c) ./ inv(c.A)[:, j]\n    epsilon[j] = minimum(eps[eps.>0])\nend\nepsilon","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We now have the smallest derease in carrying capacity that would result in the extinction of a species. The species becoming extinct can be the focal species or some other species. In the first scenario, there is no functional extinction but solely a numerical extinction.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We can compute the ecologically effective population size (EEP) of each species, that is, the minimal population size below which some species go extinct. When there is no functional extinction, but a numerical extinction, the EEP size is equal to zero.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"delta_N = Diagonal(inv(c.A)) * epsilon\nN_EEP = Neq .+ delta_N","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We see that all species EEP size are zero. This is not suprising because we have set very weak interactions strengths. Thus, the species are not dependent on each other, and the decrease in abundance of one species has very little effect on the others.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Let's repeat our analysis with a stronger interaction matrix.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"mu, sigma = -1 / S, 0.8 / sqrt(S)\nK_mean, K_std = 1, 0.3\nc = rand(\n    Community,\n    S;\n    A_ij = Normal(mu, sigma),\n    K_i = Normal(K_mean, K_std),\n    interaction = :core,\n)","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Because interactions are stronger, we expect some extinctions during the community assembly.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Neq = abundance(c)\nall(Neq .> 0)","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"So let's assemble the community, and keep only the surviving species. This can be done simply by calling the assemble function that is been designed precisely for this purpose.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"c_new = assemble(c)","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We can check that some species have gone extinct during the assembly.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"S_new = richness(c_new) # Smaller than the initial richness S.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"epsilon = zeros(S_new)\nfor j = 1:S_new\n    eps = -abundance(c_new) ./ inv(c_new.A)[:, j]\n    epsilon[j] = minimum(eps[eps.>0])\nend\ndelta_N = Diagonal(inv(c_new.A)) * epsilon\nNeq = abundance(c_new)\nN_EEP = Neq .+ delta_N","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We see that all EEP sizes are not zeros. This means that some species are functionally necessary for the communtiy.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"Finally, we can plot the EEP size of each species against their initial equilibrium abundance.","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"fig = Figure();\nax = Axis(fig[1, 1]; xlabel = \"Equilibrium abundance\", ylabel = \"EEP size\")\npositive_eep = N_EEP .> 0\nscatter!(Neq[positive_eep], N_EEP[positive_eep])\nfig","category":"page"},{"location":"functional-extinctions/","page":"Functional extinctions","title":"Functional extinctions","text":"We note a stong positive trend, that is, the more abundant the species the larger its EEP size therefore the more \"quickly\" it goes functionally extinct.","category":"page"}]
}
