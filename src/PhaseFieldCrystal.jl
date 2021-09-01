#
#  Simulate.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#

module PhaseFieldCrystal

# Import Julia packages
using DifferentialEquations
using LinearAlgebra
using Random
using NumericalIntegration
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using Plots
using JLD2

# Import local modules
include("Model.jl"); using .Model
include("CreateRunDirectory.jl"); using .CreateRunDirectory
include("InitialConditions.jl"); using .InitialConditions
include("Visualise.jl"); using .Visualise
include("ImportImage.jl"); using .ImportImage
include("FreeEnergy.jl"); using .FreeEnergy

@inline @views function phaseFieldCrystal(imagePath,L,r,ϕ₀,α₀,q,tMax,loggerFlag,outputFlag,visualiseFlag)

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

    imageMask,N = importImage(imagePath)

    # Derived parameters
    h      = L/(N-1)    # Spatial separation of grid points
    outInt = tMax/100   # Output interval
    tspan  = (0.0,tMax) # Time span for solution
    q²     = q^2        # Precalculate powers of q to reduce later calculations
    q⁴     = q^4        # Precalculate powers of q to reduce later calculations

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ = initialConditions(imageMask,L,N,α₀,ϕ₀,q)

    # Array of parameters to pass to solver
    p = [deriv, part1, part2, N, h, αᵢ, αⱼ, r, q, q², q⁴, graduᵢ, graduⱼ]

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(PFC!, u0, tspan, p)

    # Start progress logger
    loggerFlag==1 ? global_logger(TerminalLogger()) : nothing

    # Solve problem
    sol = solve(prob, alg_hints=[:stiff], reltol=10E-2, saveat=outInt, maxiters=1e9, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")

    # Calculate and plot free energy
    freeEnergies = freeEnergy(sol, N, L, q, r, h)

    if outputFlag==1
        # Create output folder and data files
        folderName = createRunDirectory(L,N,h,r,ϕ₀,α₀,q,outInt,tMax)
        # Save variables and results to file
        @info "Saving data to $folderName/data.jld2"
        jldsave("$folderName/data.jld2";sol,imageMask,freeEnergies,N,L,q,r,h,folderName)
    end

    # Plot results as animated gif
    if visualiseFlag==1 && outputFlag==1
        visualise(sol,N,h,freeEnergies,imageMask,folderName)
    end

    return 1

end

export phaseFieldCrystal

end
