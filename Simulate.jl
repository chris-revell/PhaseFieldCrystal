#
#  Simulate.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#

module Simulate

# Import Julia packages
using DifferentialEquations
using LinearAlgebra
using Random
using NumericalIntegration
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using Plots

# Import local modules
using Model
using CreateRunDirectory
using InitialConditions
using Visualise
using FreeEnergy

@inline @views function simulate(L,N,r,ϕ₀,α₀,q,tMax)

    # BLAS.set_num_threads(1)

    # Input parameters
    # L     Spatial dimensions of domain             (= 200.0 )
    # N     Number of grid points in each dimension  (= 200   )
    # r     Parameter in Swift-Hohenberg equation    (= -0.9  )
    # ϕ₀    Mean order parameter across domain       (= -0.516)
    # α     Diffusivity                              (= 1.0   )
    # q     Parameter in Swift-Hohenberg equation    (= 0.1   )
    # nPlot Number of plots produced by return       (= 100   )
    # tMax  Run time of simulation                   (= 2000.0)
    ilow  = 35
    ihigh = 70
    jlow  = 30
    jhigh = 80

    # Derived parameters
    h      = L/(N-1)    # Spatial separation of grid points
    outInt = tMax/100   # Output interval
    tspan  = (0.0,tMax) # Time span for solution
    q²     = q^2        # Precalculate powers of q to reduce later calculations
    q⁴     = q^4        # Precalculate powers of q to reduce later calculations

    # Create output folder and data files
    folderName = createRunDirectory(L,N,h,r,ϕ₀,α₀,q,outInt,tMax)

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ = initialConditions(L,N,α₀,ϕ₀,q,ilow,ihigh,jlow,jhigh)

    # Array of parameters to pass to solver
    p = [deriv, part1, part2, N, h, αᵢ, αⱼ, r, q, q², q⁴, graduᵢ, graduⱼ]

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(PFC!, u0, tspan, p)

    # Start progress logger
    global_logger(TerminalLogger())

    # Solve problem
    sol = solve(prob, alg_hints=[:stiff], reltol=10E-2, saveat=outInt, maxiters=1e9, progress=true, progress_steps=10, progress_name="PFC model")

    # Plot results as animated gif
    visualise(sol,ilow,ihigh,jlow,jhigh,N,folderName)

    freeEnergies = zeros(size(sol.u))
    for (i,u) in enumerate(sol.u)
        freeEnergies[i] = freeEnergy(u, N, L, q, r, h)
    end
    #display(plot(sol.t,freeEnergies))
    #display(freeEnergies)

    return 1

end

export simulate

end
