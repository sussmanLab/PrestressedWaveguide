@enum PropagationDirection DirX DirY
@enum Polarization Long Trans

abstract type AbstractSystemConfig end

struct LinearConfig <: AbstractSystemConfig
    kB::Float64
    S::Float64
end

struct BistableConfig <: AbstractSystemConfig
    uA::Float64
    uC::Float64
    kB::Float64 # need to use default = 1.0 for now
    S::Float64 # need to use default = 1.0 for now
end

function build_system(cfg::LinearConfig, nx::Int, ny::Int; sliding_dim::Int=0)
    return WaveguideLattice(nx, ny, cfg.kB, cfg.S; sliding_dim=sliding_dim)
end

function build_system(cfg::BistableConfig, nx::Int, ny::Int; sliding_dim::Int=0)
    return BistableWaveguideLattice(nx, ny, cfg.uA, cfg.uC; sliding_dim=sliding_dim)
end

"""
theory dispersion calculations as in the paper
"""
function theory_dispersion(::Val{DirX}, ::Val{Long}, k, kB, S)
    omega_sq = 8.0 * sin(k / 4.0)^4
    vg = sin(k / 2.0) / sqrt(2.0)
    return (omega_sq=omega_sq, vg=vg)
end

function theory_dispersion(::Val{DirX}, ::Val{Trans}, k, kB, S)
    omega_sq = 8.0 * sin(k / 4.0)^2
    vg = cos(k / 4.0) / sqrt(2.0)
    return (omega_sq=omega_sq, vg=vg)
end

function theory_dispersion(::Val{DirY}, ::Val{Long}, k, kB, S)
    omega_sq = 8.0 * sin(sqrt(3.0) * k / 4.0)^4
    vg = sqrt(3.0/2.0) * sin(sqrt(3.0) * k / 2.0)
    return (omega_sq=omega_sq, vg=vg)
end

function theory_dispersion(::Val{DirY}, ::Val{Trans}, k, kB, S)
    s_half = sqrt(3.0) * k / 2.0
    s_full = sqrt(3.0) * k
    
    omega_sq = (4.0 / 3.0) * (2.0 - 3.0 * cos(s_half) + cos(s_full))
    
    omega = omega_sq > 1e-9 ? sqrt(omega_sq) : 0.0
    
    vg_num = 3.0*sin(s_half) - 2.0*sin(s_full)
    vg_den = sqrt(3.0) * omega
    
    vg = omega > 1e-9 ? vg_num / vg_den : 0.0
    
    return (omega_sq=omega_sq, vg=vg)
end


function calc_theory(cfg::LinearConfig, dir, pol, k)
    return theory_dispersion(dir, pol, k, cfg.kB, cfg.S)
end

function calc_theory(cfg::BistableConfig, dir, pol, k)
    return theory_dispersion(dir, pol, k, cfg.kB, cfg.S) 
end

function get_param_str(cfg::LinearConfig)
    return "Harmonic_S$(cfg.S)"
end

function get_param_str(cfg::BistableConfig)
    return "Bistable_S$(cfg.S)"
end
"""
    run_dispersion_sweep(nPrimary, nSecondary, direction, polarization; kwargs...)

Runs a sweep over wavenumbers to validate dispersion relations.
Uses the WaveguideLattice constructor
"""
function run_dispersion_sweep(
    nPrimary::Int, nSecondary::Int, 
    direction::PropagationDirection, 
    polarization::Polarization,
    sys_config::AbstractSystemConfig;
    tMeasMultiples::Float64=4.0, 
    dt::Float64=1e-2, 
    amplitude::Float64=1e-4,
    error_threshold::Float64=95.0,
    output_dir::String="waveguidePropagationResults"
)
    is_x = (direction == DirX)
    is_long = (polarization == Long)
    
    case_name = (is_x ? "x" : "y") * (is_long ? "Longitudinal" : "Transverse")
    param_str = get_param_str(sys_config)
    
    filename = joinpath(output_dir, "$(case_name)_np$(nPrimary)_ns$(nSecondary)_$(param_str)_amp$(amplitude)_dt$(dt).txt")
    mkpath(output_dir)

    val_dir = is_x ? Val(DirX) : Val(DirY)
    val_pol = is_long ? Val(Long) : Val(Trans)

    open(filename, "w") do io
        writedlm(io, [["Mode", "k", "omega_theory", "vg_theory", "vg_x", "vg_y", "prho_x", "prho_y", "jE_x", "jE_y"]])

        println("\n--- $case_name Modes ---")
        @printf("%4s | %8s | %11s | %11s | %14s | %7s\n", "Mode", "k", "ω_th", "vg_th", "vg_sim", "% Err")
        println("-----|----------|-------------|-------------|----------------|--------")

        range_modes = (div(nPrimary, 2) - 1) :-1:2

        for modeIdx in range_modes
            nx = is_x ? nPrimary : nSecondary
            ny = is_x ? nSecondary : nPrimary
            
            sys, calc = build_system(sys_config, nx, ny)
            integrator = VerletIntegrator(sys, calc)

            L_dim = is_x ? sys.boundary_conditions.box_dimensions[1] : sys.boundary_conditions.box_dimensions[2]
            k = 2 * pi * modeIdx / L_dim
            
            (omega_sq, vg_theory) = calc_theory(sys_config, val_dir, val_pol, k)

            if omega_sq < 1e-9 || abs(vg_theory) < 1e-3
                continue
            end
            
            omega_theory = sqrt(omega_sq)
            period = 2 * pi / omega_theory
            tMeas = tMeasMultiples * period

            mx = is_x ? Float64(modeIdx) : 0.0
            my = is_x ? 0.0 : Float64(modeIdx)
            
            Ax, Ay = 0.0, 0.0
            if direction == DirX
                Ax = (polarization == Long) ? amplitude : 0.0
                Ay = (polarization == Trans) ? amplitude : 0.0
            else # DirY
                Ax = (polarization == Trans) ? amplitude : 0.0
                Ay = (polarization == Long) ? amplitude : 0.0
            end

            results = measure_wave_properties(
                sys, calc, integrator,
                mx, my, Ax, Ay, omega_theory,
                tMeas, dt; tSettle=0.0
            )

            measured_vg = is_x ? results.group_velocity[1] : results.group_velocity[2]
            err = 100 * abs(measured_vg - vg_theory) / abs(vg_theory)

            @printf("%4d | %8.4f | %11.4f | %11.4f | %14.4f | %6.2f%%\n", 
                    modeIdx, k, omega_theory, vg_theory, measured_vg, err)
            
            row = [modeIdx, k, omega_theory, vg_theory, 
                   results.group_velocity[1], results.group_velocity[2],
                   results.momentum_density[1], results.momentum_density[2],
                   results.energy_current_density[1], results.energy_current_density[2]]
            writedlm(io, [row])

            if err > error_threshold && modeIdx < div(nPrimary, 4) 
                println("    --> Divergence detected (>$(error_threshold)%). Halting track.")
                break
            end
        end
    end
end

"""
    measure_all_fluxes(; nPrimary=200, nSecondary=20, kwargs...)

Parallel driver to run all 4 combinations (XL, XT, YL, YT) for a set of S parameters.
"""
function measure_all_fluxes(; nPrimary=200, nSecondary=20, amplitude=1e-2, dt=1e-2)

    harmonic_configs= [LinearConfig(1.0, 1.0), LinearConfig(1.0, -1.0)]
    combinations = Iterators.product(
        [DirX, DirY], 
        [Long, Trans], 
        harmonic_configs
    )
    
    tasks = map(combinations) do (dir, pol, cfg)
        Threads.@spawn begin
            println("Starting: $dir - $pol - $(get_param_str(cfg)) on Thread $(Threads.threadid())")
            run_dispersion_sweep(
                nPrimary, nSecondary, dir, pol, cfg; 
                amplitude=amplitude, dt=dt
            )
        end
    end

    foreach(wait, tasks)
    println("All measurements complete.")
end

function measure_bistable_fluxes(; nPrimary=200, nSecondary=20, amplitude=1e-2, dt=1e-2)

    # uA and uC chosen so that δ ~ 0.75 of r0?
    bistable_configs= [BistableConfig(0.0703125,0.7,1.0, 1.0)]
    combinations = Iterators.product(
        [DirX, DirY], 
        [Long, Trans], 
        bistable_configs
    )
    
    tasks = map(combinations) do (dir, pol, cfg)
        Threads.@spawn begin
            println("Starting Bistable: $dir - $pol - $(get_param_str(cfg)) on Thread $(Threads.threadid())")
            run_dispersion_sweep(
                nPrimary, nSecondary, dir, pol, cfg; 
                amplitude=amplitude, dt=dt
            )
        end
    end

    foreach(wait, tasks)
    println("All measurements complete.")
end

"""
    run_sliding_dispersion_sweep(nPrimary, nSecondary, direction, sys_config; kwargs...)

Runs a sweep over wavenumbers for longitudinal waves with sliding boundary conditions.
"""
function run_sliding_dispersion_sweep(
    nPrimary::Int, nSecondary::Int, 
    direction::PropagationDirection, 
    sys_config::AbstractSystemConfig;
    tMeasMultiples::Float64=4.0, 
    dt::Float64=1e-2, 
    amplitude::Float64=1e-4,
    error_threshold::Float64=95.0,
    output_dir::String="slidingWaveguidePropagationResults"
)
    is_x = (direction == DirX)
    polarization = Long
    
    case_name = (is_x ? "x" : "y") * "Longitudinal_Sliding"
    param_str = get_param_str(sys_config)
    
    filename = joinpath(output_dir, "$(case_name)_np$(nPrimary)_ns$(nSecondary)_$(param_str)_amp$(amplitude)_dt$(dt).txt")
    mkpath(output_dir)

    val_dir = is_x ? Val(DirX) : Val(DirY)
    val_pol = Val(Long)
    
    # sliding_dim is orthogonal to direction
    s_dim = is_x ? 2 : 1

    open(filename, "w") do io
        writedlm(io, [["Mode", "k", "omega_theory", "vg_theory", "vg_x", "vg_y", "prho_x", "prho_y", "jE_x", "jE_y"]])

        println("\n--- $case_name Modes ---")
        @printf("%4s | %8s | %11s | %11s | %14s | %7s\n", "Mode", "k", "ω_th", "vg_th", "vg_sim", "% Err")
        println("-----|----------|-------------|-------------|----------------|--------")

        range_modes = (div(nPrimary, 2) - 1) :-1:2

        for modeIdx in range_modes
            nx = is_x ? nPrimary : nSecondary
            ny = is_x ? nSecondary : nPrimary
            
            sys, calc = build_system(sys_config, nx, ny; sliding_dim=s_dim)
            integrator = VerletIntegrator(sys, calc)

            L_dim = is_x ? sys.boundary_conditions.box_dimensions[1] : sys.boundary_conditions.box_dimensions[2]
            k = 2 * pi * modeIdx / L_dim
            
            (omega_sq, vg_theory) = calc_theory(sys_config, val_dir, val_pol, k)

            if omega_sq < 1e-9 || abs(vg_theory) < 1e-3
                continue
            end
            
            omega_theory = sqrt(omega_sq)
            period = 2 * pi / omega_theory
            tMeas = tMeasMultiples * period

            mx = is_x ? Float64(modeIdx) : 0.0
            my = is_x ? 0.0 : Float64(modeIdx)
            
            Ax, Ay = 0.0, 0.0
            if direction == DirX
                Ax = amplitude
                Ay = 0.0
            else
                Ax = 0.0
                Ay = amplitude
            end

            results = measure_wave_properties(
                sys, calc, integrator,
                mx, my, Ax, Ay, omega_theory,
                tMeas, dt; tSettle=0.0
            )

            measured_vg = is_x ? results.group_velocity[1] : results.group_velocity[2]
            err = 100 * abs(measured_vg - vg_theory) / abs(vg_theory)

            @printf("%4d | %8.4f | %11.4f | %11.4f | %14.4f | %6.2f%%\n", 
                    modeIdx, k, omega_theory, vg_theory, measured_vg, err)
            
            row = [modeIdx, k, omega_theory, vg_theory, 
                   results.group_velocity[1], results.group_velocity[2],
                   results.momentum_density[1], results.momentum_density[2],
                   results.energy_current_density[1], results.energy_current_density[2]]
            writedlm(io, [row])

            if err > error_threshold && modeIdx < div(nPrimary, 4) 
                println("    --> Divergence detected (>$(error_threshold)%). Halting track.")
                break
            end
        end
    end
end

"""
    measure_sliding_fluxes(; nPrimary=200, nSecondary=20, kwargs...)

Parallel driver to run longitudinal combinations (XL, YL) with sliding boundaries.
"""
function measure_sliding_fluxes(; nPrimary=200, nSecondary=20, amplitude=1e-2, dt=1e-2)
    harmonic_configs = [LinearConfig(1.0, 1.0), LinearConfig(1.0, -1.0)]
    combinations = Iterators.product(
        [DirX, DirY], 
        harmonic_configs
    )
    
    tasks = map(combinations) do (dir, cfg)
        Threads.@spawn begin
            println("Starting Sliding: $dir - Long - $(get_param_str(cfg)) on Thread $(Threads.threadid())")
            run_sliding_dispersion_sweep(
                nPrimary, nSecondary, dir, cfg; 
                amplitude=amplitude, dt=dt
            )
        end
    end

    foreach(wait, tasks)
    println("All sliding measurements complete.")
end

