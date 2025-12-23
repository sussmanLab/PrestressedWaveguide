"""
    run_simulation!(system, force_calculator, integrator, dt, num_steps; callback, every_n_steps=1)

Runs a complete N-body simulation for a specified number of discrete time steps.

This is the main simulation loop. It repeatedly calls the `integrate!` method,
which dispatches on the provided `integrator` type to advance the `system` state.
The `system` object is mutated in-place over the course of the simulation.
Meanwhile, a `callback` function whose signature must be `f(system,step_number)` is called every `every_n_steps`.
By default, a no-op function is used as the callback at every step

# Arguments
- `system::System{D,T}`: The system of particles to be simulated. This object is mutated.
- `force_calculator::AbstractForceCalculator`: A ForceCalculator object that holds an algorithm (e.g., O(N^2) brute force) and a force law that calculates the force between two particles. 
- `integrator::AbstractIntegrator`: An integrator object (e.g., `VerletIntegrator`) that holds the algorithm and its state. *Note:* If this object requires initialization, that must be done before being passed to this function.
- `dt::Float64`: The size of a single time step.
- `num_steps::Int`: The total number of simulation steps to perform.

# Keyword Arguments
- `callback`: A function that is called every `every_n_steps`. It should have the signature `f(system, step_number)`.
- `every_n_steps::Int=1`: The frequency at which to call the callback function.

# Returns
- `nothing`: This function modifies the `system` object in-place.
"""
function run_simulation!(system,force_calculator,integrator,time_step_size,number_of_steps; callback = (sys, step) -> nothing, every_n_steps::Int=1)

    for i in 1:number_of_steps
        integrate!(system,force_calculator,integrator,time_step_size)
        if i % every_n_steps == 0
            callback(system, i)
        end
    end
end

