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
include("CreateLaplacian.jl"); using .CreateLaplacian
include("CreateRunDirectory.jl"); using .CreateRunDirectory
include("InitialConditions.jl"); using .InitialConditions
include("Visualise.jl"); using .Visualise
include("FreeEnergy.jl"); using .FreeEnergy

function phaseFieldCrystal(nGrid,lSpace,r,ϕ₀,a,dt,tMax,loggerFlag,outputFlag,visualiseFlag,integrator,randomOrNot)

    #BLAS.set_num_threads(8)

    # Input parameters
    # nGrid         Number of grid points in both dimensions  (eg. = 200          )
    # lSpace        Spatial dimensions of domain              (eg. = 200.0        )
    # r             Parameter in Swift-Hohenberg equation     (eg. = -0.9 or 0.5  )
    # ϕ₀            Mean order parameter across domain        (eg. = -0.516       )
    # a             Parameter in splitting scheme             (eg. = 2.0          )
    # dt            Time step for implicit semi-linear scheme (eg. = 0.01         )
    # tMax          Run time of simulation                    (eg. = 20.0         )
    # loggerFlag    Display progress bar in REPL when =1; Always set to 0 on HPC
    # outputFlag    Flag to control whether data are saved to file (=1 or 0)
    # visualiseFlag Flag to control whether data are plotted  (=1 or 0)
    # integrator    Controls which integration scheme to use ="split" or "explicit"

    
    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,mat1,mat2,h = initialConditions(lSpace,nGrid,ϕ₀,1.0,randomOrNot)

    ∇² = createLaplacian(nGrid,h)
    
    linearOperator = (1.0-r+a).*∇² .+ ∇²*∇²*∇²

    # Array of parameters to pass to solver
    p = [∇², linearOperator, mat1, mat2, r, a]

    # Start progress logger if loggerFlag argument is 1
    loggerFlag==1 ? global_logger(TerminalLogger()) : nothing

    # Define ODE problem using discretised derivatives
    if integrator=="split"
        prob = SplitODEProblem(DiffEqArrayOperator(linearOperator),splitNonlinearPart!,u0,(0.0,tMax),p)
        sol = solve(prob, LawsonEuler(krylov=true, m=50), dt=0.01, saveat=(tMax/100), rel_tol=0.001, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")
    elseif integrator=="explicit"
        prob = ODEProblem(PFC!,u0,(0.0,tMax),p)
        sol = solve(prob, Tsit5(), saveat=(tMax/100), rel_tol=0.001, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")
    end 
    #prob = SplitODEProblem(splitLinearPart!,splitNonlinearPart!,u0,(0.0,tMax),p)
    #sol = solve(prob, IMEXEuler(), dt=0.0001, saveat=(tMax/100), rel_tol=0.001, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")
    
    # Calculate and plot free energy
    freeEnergies = freeEnergy(sol, ∇², mat1, mat2, nGrid, lSpace, r)

    if outputFlag==1
        # Create output folder and data files
        folderName = createRunDirectory(lSpace,nGrid,h,r,ϕ₀,tMax/100,tMax)
        # Save variables and results to file
        @info "Saving data to $folderName/data.jld2"
        jldsave("$folderName/data.jld2";sol,∇²,freeEnergies,nGrid,lSpace,r,h,folderName)
        # Plot results as animated gif
        if visualiseFlag==1
            visualise(sol,∇²,nGrid,freeEnergies,folderName)
        end
    end

    return 1

end

export phaseFieldCrystal

end
