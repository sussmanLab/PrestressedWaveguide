module PrestressedWaveguide

using StaticArrays, Random, Printf
using Base.Threads
using DelimitedFiles

export norm_squared, spring_force, HarmonicSpring, BistableQuarticSpring, Bond, potential_energy
export Particle, ParticleSystem, PeriodicBoundaryConditions, SlidingWaveguideBCs, LongitudinalConstraintBC, box_displacement, apply_all_bcs!, kinetic_energy, total_energy
export VerletIntegrator, FixedTopologyCalculator, compute_acceleration!, integrate!
export OneDChainLattice, SimpleTriangularLattice, WaveguideLattice, BistableWaveguideLattice
export imposeTravelingWavePacket!, imposeTravelingWave!, imposeWave!
export run_simulation!, measure_wave_properties, calculate_pseudomomentum_density, average_energy_density
export DirX, Long, Trans, DirY, LinearConfig, BistableConfig
export run_dispersion_sweep, measure_all_fluxes, measure_bistable_fluxes, run_sliding_dispersion_sweep, measure_sliding_fluxes, run_longitudinal_constraint_dispersion_sweep, measure_longitudinal_constraint_fluxes

include("forceLaws.jl")

include("particleSystemStructs.jl")
include("forceCalculators.jl")

include("latticeConstructors.jl")

include("integrators.jl")
include("simulator.jl")

include("waveMeasurements.jl")

include("prestressedDispersionSweeps.jl")

end # module PrestressedWaveguide
