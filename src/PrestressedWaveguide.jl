module PrestressedWaveguide

using StaticArrays, Random
using Base.Threads
using DelimitedFiles

export norm_squared, spring_force, HarmonicSpring, BistableQuarticSpring, Bond, potential_energy
export Particle, ParticleSystem, PeriodicBoundaryConditions, box_displacement, apply_all_bcs!, kinetic_energy, total_energy
export VerletIntegrator, FixedTopologyCalculator, compute_acceleration!, integrate!
export OneDChainLattice, SimpleTriangularLattice, WaveguideLattice, imposeTravelingWavePacket!, imposeTravelingWave!, imposeWave!
export run_simulation!

include("forceLaws.jl")

include("particleSystemStructs.jl")
include("forceCalculators.jl")

include("latticeConstructors.jl")

include("integrators.jl")
include("simulator.jl")

end # module PrestressedWaveguide
