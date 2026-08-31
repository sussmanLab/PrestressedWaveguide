#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using PrestressedWaveguide
Pkg.activate("plotting"; shared=true)
using GLMakie
using StaticArrays
using LinearAlgebra
using Printf

const NX = 60           
const NY = 20          
const KB = 1.0          
const SXX = 0.0         
const DT = 0.01      
const FRAMES = 200     
const STEPS_PER_FRAME = 500 

const MODE_q = 3.0      
const AMP = 0.2        

"""
Helper function to get bond segments for plotting.
Crucially, this skips bonds that wrap around Periodic Boundary Conditions
so we don't get ugly lines drawing across the entire screen.
"""
function get_bond_segments(system, bonds, box_size)
    segments = Point2f[]
    Lx, Ly = box_size[1], box_size[2]
    
    for bond in bonds
        p1 = system.particles[bond.i].position
        p2 = system.particles[bond.j].position
        
        dx = abs(p1[1] - p2[1])
        dy = abs(p1[2] - p2[2])
        
        if dx < Lx/2 && dy < Ly/2
            push!(segments, Point2f(p1[1], p1[2]))
            push!(segments, Point2f(p2[1], p2[2]))
        end
    end
    return segments
end

function animateSpinWave_experiment()
    println("Initializing System...")
    
    sys_config= PrestressedWaveguide.LinearConfig(KB, SXX)
    sys, calc = PrestressedWaveguide.build_system(sys_config, NX, NY)

    L_dims = sys.boundary_conditions.box_dimensions
    k = 2 * pi * MODE_q / L_dims[1]
    
    # from Eq. S52,53 when S=0
    s_half = sin(k / 2.0)
    c_half = cos(k / 2.0)
    omega_sq = 2.0 * KB * ( -s_half^2 + 2.0 * (1.0 - c_half) )
    omega = sqrt(omega_sq)
    
    println("Imposing Spin Wave (k=$k, ω=$omega)...")
    
    imposeSpinWave!(sys, MODE_q, AMP, omega)
    integrator = VerletIntegrator(sys, calc)

    println("Setting up GLMakie...")
    
    points_obs = Observable([Point2f(p.position...) for p in sys.particles])
    
    colors_obs = Observable([p.velocity[2] for p in sys.particles])
    
    bonds_obs = Observable(get_bond_segments(sys, calc.bonds, L_dims))

    fig = Figure(size = (1000, 400), fontsize=18)
    ax = Axis(fig[1, 1], 
              title = "Spin Wave Animation", 
              xlabel = "x", ylabel = "y",
              aspect = DataAspect())
    
    xlims!(ax, -1, L_dims[1] + 1)
    ylims!(ax, -1, L_dims[2] + 1)

    linesegments!(ax, bonds_obs, color = :gray80, linewidth = 1.0)

    sc = scatter!(ax, points_obs, 
             color = colors_obs, 
             colormap = :RdBu, 
             colorrange = (-AMP * omega, AMP * omega), # the other script should probably be this, too -- map colors to velocity instead of displacement
             markersize = 8)
    Colorbar(fig[1, 2], sc, label = "vertical velocity (vy)")

    println("Starting Recording...")
    output_filename = "spin_wave_S$(SXX).mp4"
    
    record(fig, output_filename, 1:FRAMES; framerate = 30) do frame
        run_simulation!(sys, calc, integrator, DT, STEPS_PER_FRAME)
        
        points_obs[] = [Point2f(p.position...) for p in sys.particles]
        colors_obs[] = [p.velocity[2] for p in sys.particles]
        bonds_obs[] = get_bond_segments(sys, calc.bonds, L_dims)
        
        if frame % 50 == 0
            print("Rendered frame $frame / $FRAMES \r")
        end
    end

end

animateSpinWave_experiment()
