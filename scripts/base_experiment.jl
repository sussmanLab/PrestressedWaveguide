#!/usr/bin/env julia

using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

Pkg.instantiate()

using PrestressedWaveguide
using Dates

const N_PRIMARY = 250
const N_SECONDARY = 20
const AMPLITUDE = 1e-4
const DT1 = 1e-2
const DT2 = 1e-4

println("==================================================")
println("Experiment Started at: ", now())
println("Parameters:")
println("  nPrimary   : ", N_PRIMARY)
println("  nSecondary : ", N_SECONDARY)
println("  Amplitude  : ", AMPLITUDE)
println("  Time Step  : ", DT1, " ", DT2)
println("==================================================")

try
    measure_all_fluxes(; 
        nPrimary = N_PRIMARY, 
        nSecondary = N_SECONDARY,
        amplitude = AMPLITUDE,
        dt = DT1
    )
    measure_all_fluxes(; 
        nPrimary = N_PRIMARY, 
        nSecondary = N_SECONDARY,
        amplitude = AMPLITUDE,
        dt = DT2
    )
    println("\nSUCCESS: All experiments completed successfully.")
catch e
    println("\nERROR: Experiment failed with error:")
    showerror(stdout, e)
    exit(1)
end

println("Finished at: ", now())
