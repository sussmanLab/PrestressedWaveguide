# PrestressedWaveguide.jl

This repo contains the simulation and analysis code in support of "Conservative yet constitutively odd elasticity in prestressed metamaterials", Tyler A. Engstrom and Daniel M. Sussman. Accepted for publication in PRL as of Nov 16, 2026 (preprint available [on the arXiv:2509.19560](https://arxiv.org/abs/2509.19560)).

## Paper overview and findings

Our paper looks at nonlinear continuum mechanics of solids perturbed around an equilibrium prestressed reference state. We examine the effective constitutive tensor $B_{ijkl}$ that relates variations in the Cauchy stress $\sigma_{ij}$ to the Lagrangian strain $\epsilon_{kl}$:

$$\sigma_{ij} = S_{ij} + B_{ijkl}\epsilon_{kl} + \text{rotation-dependent terms}$$

Because $\sigma_{ij}$ and $\epsilon_{kl}$ do not form an energy-conjugate pair, $B_{ijkl}$ does not need to possess major symmetry even when the underlying potential energy is fully conservative. Its anti-symmetric part is governed directly by the equilibrium Cauchy prestress $S_{ij}$:

$$B_{ijkl} - B_{klij} = S_{kl}\delta_{ij} - S_{ij}\delta_{kl}$$

By applying a specific anisotropic prestress ($S_{xx} = -S_{yy}$), we inverse-design 2D Bravais lattice metamaterials in which the symmetric part of $B_{ijkl}$ vanishes entirely, driving the system to the purely anti-symmetric limit $B_{ijkl} = -B_{klij}$. Does this limit correspond to anything physically interesting about the system? Our papers explores the answer to that question (demonstrating things like "DC" energy transport, chiral phonons, and in-plane longitudinal soft modes with a quadratic dispersion.

## Simulations and repository structure

This Julia pacakge implements an MD engine in the NVE ensemble to test the continuum predictions for finite-amplitude waves.

### Core Package (`src/`)

Some of the core files and the things they implement:

* Metamaterial Lattice Construction (`latticeConstructors.jl`, `forceLaws.jl`): Utilities for constructing 2D triangular lattice waveguides with nearest- and next-nearest-neighbor interactions partitioned into three families of different kinds of springs (see the paper for details). We have functions for studying both `HarmonicSpring` systems and springs made out of `BistableQuarticSpring` elements matched to the target equilibrium forces and curvatures.

* Symplectic Integration and Boundary Conditions (`integrators.jl`, `particleSystemStructs.jl`, `simulator.jl`): Contains standard code for evolving particle trajectories using the Velocity Verlet algorithm. To test the sensitivity of long-wavelength modes to off-axis instabilities, we have looked at a variety of boundary conditions: `PeriodicBoundaryConditions`, `SlidingWaveguideBCs`, and `LongitudinalConstraintBC`.

* Transport and Wave Measurements (`waveMeasurements.jl`): Functions for computing instantaneous and time-averaged energy fluxes across fixed planes (using the standard Irving-Kirkwood formulation), evaluating lattice pseudomomentum density from local deformation gradients, and extracting simulated group velocities.

* Automated Dispersion Sweeps (`prestressedDispersionSweeps.jl`): Functions that orchestrate some of the sweeps across lattice wavenumbers, $q$ and evelove the system over multiple wave periods.


## How to Generate Data

To reproduce the primary molecular dynamics dataset from the paper (sweeping over wavevectors $q$, propagation directions, polarizations, and prestress values $S$), clone this repository and run `scripts/varyStressSweeps.jl` from the repository root (optionally, with with multi-threading enabled):

```shell
julia --threads auto scripts/varyStressSweeps.jl
```

This script automatically activates and instantiates the project environment and writes the output tables to the `waveguidePropagationResults/` directory.

## Additional experiments and visualizations

The `scripts/` directory contains a few additional scripts that may be of interest. 

### Animations

Two scripts will render `.mp4` videos of traveling waves, or circularly polarized spin waves at $S = 0$:
```shell
julia scripts/animation_experiment.jl
julia scripts/spinWave_experiment.jl
```
These make use of DMS' shared `@plotting` environment -- to generate, you should add `GLMakie` to a shared environment of this name on your installation.

# Experiments

To test mode stability when negative-stiffness harmonic springs are replaced by quartic double-well potentials:
  ```shell
  julia --threads auto scripts/bistable_experiment.jl
  julia --threads auto scripts/varyStressSweeps.jl
  ```
To run longitudinal dispersion sweeps under sliding waveguide or longitudinally constrained boundary conditions:
  ```shell
  julia --threads auto scripts/slidingBC_experiment.jl
  julia --threads auto scripts/longitudinalConstraintBC_experiment.jl
  ```
