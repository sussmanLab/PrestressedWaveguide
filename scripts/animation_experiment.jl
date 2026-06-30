#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using PrestressedWaveguide
Pkg.activate("plotting"; shared=true)
using GLMakie
using StaticArrays
using LinearAlgebra
using Printf

const NX = 20           
const NY = 200          
const KB = 1.0          
const SXX = 1.0         
const DT = 0.001      
const FRAMES = 300     
const STEPS_PER_FRAME = 1000 

const MODE_q = 10.0      
const AMP = 0.01        
const SIGMA = 8.0      

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

function animate_experiment()
    println("Initializing System...")
    
    sys_config= PrestressedWaveguide.BistableConfig(0.0703125,0.7,1.0, 1.0)
    sys, calc = PrestressedWaveguide.build_system(sys_config, NX, NY)
    integrator = VerletIntegrator(sys, calc)
    dir = Val(PrestressedWaveguide.DirY)
    pol = Val(PrestressedWaveguide.Long)

    L_dims = sys.boundary_conditions.box_dimensions
    k = 2 * pi * MODE_q / L_dims[2]
    (omega_sq, vg) = PrestressedWaveguide.calc_theory(sys_config,dir, pol, k)
    omega = sqrt(omega_sq)
    
    println("Imposing Wave Packet (k=$k, ω=$omega)...")
    
    imposeTravelingWave!(sys, 0.0, MODE_q, 0.0, AMP, omega)

    println("Setting up GLMakie...")
    
    points_obs = Observable([Point2f(p.position...) for p in sys.particles])
    
    colors_obs = Observable([p.velocity[2] for p in sys.particles])
    
    bonds_obs = Observable(get_bond_segments(sys, calc.bonds, L_dims))

    fig = Figure(size = (400, 1000), fontsize=18)
    ax = Axis(fig[1, 1], 
              title = "Prestressed Waveguide Animation", 
              xlabel = "x", ylabel = "y",
              aspect = DataAspect())
    
    xlims!(ax, -1, L_dims[1] + 1)
    ylims!(ax, -1, L_dims[2] + 1)

    linesegments!(ax, bonds_obs, color = :gray80, linewidth = 1.0)

    sc = scatter!(ax, points_obs, 
             color = colors_obs, 
             colormap = :RdBu, 
             colorrange = (-AMP, AMP),
             markersize = 8)

    Colorbar(fig[1, 2], sc, label = "vertical velocity (vy)")

    println("Starting Recording...")
    output_filename = "waveguide_packet_S$(SXX).mp4"
    
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

animate_experiment()
