# GLV.jl

A package for simulating the dynamics of the Generalized Lotka-Volterra (GLV) model in Julia.

## Features

- Efficient numerical integration using Julia’s differential equation solvers.
- Implement predictions for community properties using the cavity method.

For detailed usage and examples, please refer to the [documentation](https://ismael-lajaaiti.github.io/GLV.jl/).

## TODO

- [ ] Implement species sensitivity to press
- [ ] Implement species sensitivity to extinction
- [ ] Implement jacobian, reactivity, resilience, etc.
- [ ] Document stability analysis at the community
- [ ] Add a dispersion term in the model.
- [ ] Publish the package in the Julia registry.
- [x] Document species sensitivity to pulse (controlled by relative yield)
- [x] Check and document implementation of the cavity predictions
- [x] Implement functions from the cavity method to predict community properties.
- [x] Add perturbations functions (noise, extinction, press, pulse).
- [x] Implement species reactivity

