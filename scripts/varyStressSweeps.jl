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
const S_VALUES = [1.0, 0.3, 0.1, 0.03, 0.01, 0.001, 0.0001]

"""
    compute_bistable_u(L, K, F; target_ratio=0.9)

Calculates the required bistable potential depth `u` to match a specific 
length (L), stiffness (K), and force (F), while guaranteeing that the 
separation delta remains a safe fraction of the central distance r0 
(defaulting to 75%).
"""
function compute_bistable_u(L, K, F; target_ratio=0.9)
    if abs(F) < 1e-12
        return (target_ratio * L)^2 * (-K) / 4.0
    end

    R2 = target_ratio^2
    
    c3 = K * (1.0 - R2)
    c2 = -3.0*F + R2*(F + 2.0*L*K)
    c1 = -R2 * L * (L*K + 2.0*F)
    c0 = R2 * L^2 * F
    
    g(x) = c3*x^3 + c2*x^2 + c1*x + c0
    
    x_left = min(0.0, F/K)
    x_right = max(0.0, F/K)
    v_left = g(x_left)
    
    for _ in 1:60
        x_mid = (x_left + x_right) / 2.0
        v_mid = g(x_mid)
        if sign(v_mid) == sign(v_left)
            x_left = x_mid
            v_left = v_mid
        else
            x_right = x_mid
        end
    end
    
    x = (x_left + x_right) / 2.0
    
    u = (x * (K*x - 3.0*F)^2) / (8.0 * (K*x - F))
    return u
end

"""
    get_bistable_config(S; kB=1.0)

Generates a perfectly scaled BistableConfig for a given prestress S.
"""
function get_bistable_config(S; kB=1.0)
    KA = -kB / 2.0
    KC = -kB / 2.0
    
    LA = 1.0
    LC = sqrt(3.0)
    
    rA_harm = S
    rC_harm = -S / sqrt(3.0)
    
    FA = KA * (LA - rA_harm)
    FC = KC * (LC - rC_harm)
    
    uA = compute_bistable_u(LA, KA, FA)
    uC = compute_bistable_u(LC, KC, FC)
    
    return BistableConfig(uA, uC, kB, S)
end

"""
    measure_qx_fluxes(; kwargs...)

Custom parallel driver to run only the DirX (qx) combinations for both 
Longitudinal and Transverse polarizations, for both Harmonic and Bistable lattices.
"""
function measure_qx_fluxes(; nPrimary=200, nSecondary=20, amplitude=1e-2, dt=1e-3, S_values=[1.0, -1.0])
    
    harmonic_configs = [LinearConfig(1.0, S) for S in S_values]
    
    bistable_configs = [get_bistable_config(S) for S in S_values]

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
