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

@inline @views function phaseFieldCrystal(N,L,r,ϕ₀,a,tMax,loggerFlag,outputFlag,visualiseFlag)

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

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,mat1,mat2,mat3,h = initialConditions(L,N,ϕ₀,1.0)

    linearOperator,laplacianMatrix = f1(N,r,a,h)

    # Array of parameters to pass to solver
    p = [laplacianMatrix, mat1, mat2, mat3, N, h, r, a]

    # Define ODE problem using discretised derivatives
    prob = SplitODEProblem(linearOperator,f2!,u0,(0.0,tMax),p)

    # Start progress logger
    loggerFlag==1 ? global_logger(TerminalLogger()) : nothing

    # Solve problem
    sol = solve(prob, NorsettEuler(), dt=(tMax/10000.0), saveat=(tMax/100), rel_tol=0.001, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")

    # Calculate and plot free energy
    freeEnergies = freeEnergy(sol, laplacianMatrix, mat1, mat2, N, L, r)

    if outputFlag==1
        # Create output folder and data files
        folderName = createRunDirectory(L,N,h,r,ϕ₀,tMax/100,tMax)
        # Save variables and results to file
        @info "Saving data to $folderName/data.jld2"
        jldsave("$folderName/data.jld2";sol,freeEnergies,N,L,r,h,folderName)
    end

    # Plot results as animated gif
    if visualiseFlag==1 && outputFlag==1
        visualise(sol,laplacianMatrix,mat1,N,h,freeEnergies,folderName)
    end

    return 1

end

export phaseFieldCrystal

end
