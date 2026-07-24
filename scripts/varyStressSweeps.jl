#!/usr/bin/env julia

using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PrestressedWaveguide
using PrestressedWaveguide: run_dispersion_sweep, LinearConfig, BistableConfig, DirX, Long, Trans
using Dates

const N_PRIMARY = 200
const N_SECONDARY = 20
const AMPLITUDE = 1e-2
const DT = 1e-3
const S_VALUES = [1.0, 0.1, 0.01]

"""
    measure_qx_fluxes(; kwargs...)

Custom parallel driver to run only the DirX (qx) combinations for both 
Longitudinal and Transverse polarizations, for both Harmonic and Bistable lattices.
"""
function measure_qx_fluxes(; nPrimary=200, nSecondary=20, amplitude=1e-2, dt=1e-3, S_values=[1.0, -1.0])
    
    harmonic_configs = [LinearConfig(1.0, S) for S in S_values]
    
    bistable_configs = [BistableConfig(0.0703125, 0.7, 1.0, S) for S in S_values]

    configs_to_run = vcat(harmonic_configs, bistable_configs)
    polarizations = [Long, Trans]

    combinations = Iterators.product(polarizations, configs_to_run)
    
    tasks = map(combinations) do (pol, cfg)
        Threads.@spawn begin
            sys_type = cfg isa LinearConfig ? "Harmonic" : "Bistable"
            println("Starting: DirX - $pol - $sys_type (S=$(cfg.S)) on Thread $(Threads.threadid())")
            
            run_dispersion_sweep(
                nPrimary, nSecondary, DirX, pol, cfg; 
                amplitude=amplitude, dt=dt
            )
        end
    end

    foreach(wait, tasks)
    println("All DirX (qx) measurements complete.")
end

println("==================================================")
println("Experiment Started at: ", now())
println("Parameters:")
println("  nPrimary   : ", N_PRIMARY)
println("  nSecondary : ", N_SECONDARY)
println("  Amplitude  : ", AMPLITUDE)
println("  Time Step  : ", DT)
println("  S Values   : ", S_VALUES)
println("==================================================")

try
    measure_qx_fluxes(; 
        nPrimary = N_PRIMARY, 
        nSecondary = N_SECONDARY,
        amplitude = AMPLITUDE,
        dt = DT,
        S_values = S_VALUES
    )
    println("\nSUCCESS: All experiments completed successfully.")
catch e
    println("\nERROR: Experiment failed with error:")
    showerror(stdout, e)
    exit(1)
end

println("Finished at: ", now())
