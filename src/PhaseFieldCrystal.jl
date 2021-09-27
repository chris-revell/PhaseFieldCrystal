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
using SparseArrays

# Import local modules
include("Model.jl"); using .Model
include("Laplacian.jl"); using .Laplacian
include("CreateRunDirectory.jl"); using .CreateRunDirectory
include("InitialConditions.jl"); using .InitialConditions
include("Visualise.jl"); using .Visualise
include("FreeEnergy.jl"); using .FreeEnergy

@inline @views function phaseFieldCrystal(nGrid,lSpace,r,ϕ₀,a,tMax,loggerFlag,outputFlag,visualiseFlag)

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
    nGhosts = 6
    
    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,mat1,mat2,mat3,h = initialConditions(lSpace,nGrid,nGhosts,ϕ₀,1.0)

    ∇² = createLaplacian(nGrid+nGhosts,h)

    linearOperator = DiffEqArrayOperator((1.0-r+a).*∇² .+ ∇²*∇²*∇²)

    # Array of parameters to pass to solver
    p = [∇², mat1, mat2, mat3, nGrid, h, r, a]

    # Define ODE problem using discretised derivatives
    prob = SplitODEProblem(linearOperator,f2!,u0,(0.0,tMax),p)

    # Start progress logger if loggerFlag argument is 1
    loggerFlag==1 ? global_logger(TerminalLogger()) : nothing

    # Solve problem
    sol = solve(prob, SplitEuler(), dt=(tMax/100000.0), saveat=(tMax/100), rel_tol=0.001, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")

    # Calculate and plot free energy
    freeEnergies = freeEnergy(sol, ∇², mat1, mat2, nGrid, nGhosts, lSpace, r)

    if outputFlag==1
        # Create output folder and data files
        folderName = createRunDirectory(lSpace,nGrid,h,r,ϕ₀,tMax/100,tMax)
        # Save variables and results to file
        @info "Saving data to $folderName/data.jld2"
        jldsave("$folderName/data.jld2";sol,freeEnergies,nGrid,nGhosts,lSpace,r,h,folderName)
    end

    # Plot results as animated gif
    if visualiseFlag==1 && outputFlag==1
        visualise(sol,∇²,mat1,nGrid,nGhosts,h,freeEnergies,folderName)
    end

    return 1

end

export phaseFieldCrystal

end
