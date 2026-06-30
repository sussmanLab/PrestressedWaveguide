#!/usr/bin/env julia

using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

Pkg.instantiate()

using PrestressedWaveguide
using Dates

const N_PRIMARY = 200
const N_SECONDARY = 20
const AMPLITUDE = 1e-2
const DT1 = 1e-3
const DT2 = 1e-4
const S_VALUES = [-1.0, 1.0, 0.1, 0.01, -0.1, -0.01]

println("==================================================")
println("Experiment Started at: ", now())
println("Parameters:")
println("  nPrimary   : ", N_PRIMARY)
println("  nSecondary : ", N_SECONDARY)
println("  Amplitude  : ", AMPLITUDE)
println("  Time Step  : ", DT1, " ", DT2)
println("  S Values   : ", S_VALUES)
println("==================================================")

try
    measure_all_fluxes(; 
        nPrimary = N_PRIMARY, 
        nSecondary = N_SECONDARY,
        amplitude = AMPLITUDE,
        dt = DT1,
        S_values = S_VALUES
    )
    measure_all_fluxes(; 
        nPrimary = N_PRIMARY, 
        nSecondary = N_SECONDARY,
        amplitude = AMPLITUDE,
        dt = DT2,
        S_values = S_VALUES
    )
    println("\nSUCCESS: All experiments completed successfully.")
catch e
    println("\nERROR: Experiment failed with error:")
    showerror(stdout, e)
    exit(1)
end

println("Finished at: ", now())
