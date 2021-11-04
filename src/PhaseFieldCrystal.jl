#
#  PhaseFieldCrystal.jl
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
using DrWatson

# Import local modules
include("Model.jl"); using .Model
include("CreateLaplacian.jl"); using .CreateLaplacian
include("CreateGrad.jl"); using .CreateGrad
include("CreateDiv.jl"); using .CreateDiv
include("InitialConditions.jl"); using .InitialConditions
include("Visualise.jl"); using .Visualise
include("FreeEnergy.jl"); using .FreeEnergy

function phaseFieldCrystal(nGrid,lSpace,r,ϕ0,a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)

    BLAS.set_num_threads(nBlasThreads)

    # Input parameters
    # nGrid         Number of grid points in both dimensions  (eg. = 200          )
    # lSpace        Spatial dimensions of domain              (eg. = 200.0        )
    # r             Parameter in Swift-Hohenberg equation     (eg. = -0.9 or 0.5  )
    # ϕ0            Mean order parameter across domain        (eg. = -0.516       )
    # a             Parameter in splitting scheme             (eg. = 2.0          )
    # δt            Time step for implicit semi-linear scheme (eg. = 0.01         )
    # tMax          Run time of simulation                    (eg. = 20.0         )
    # loggerFlag    Display progress bar in REPL when =1; Always set to 0 on HPC
    # outputFlag    Flag to control whether data are saved to file (=1 or 0)
    # visualiseFlag Flag to control whether data are plotted  (=1 or 0)
    # integrator    Controls which integration scheme to use ="split" or "explicit"

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,mat1,mat2,h,α   = initialConditions(lSpace,nGrid,ϕ0,1.0,1)

    # Create finite difference matrices for given system parameters
    ∇² = createLaplacian(nGrid,h)
    divalphagrad = createGrad(nGrid,h,α)

    # Create matrix for linaer component of PFC equation
    linearOperator = divalphagrad.*(1.0-r+a) .+ divalphagrad*∇²*∇²

    # Array of parameters to pass to solver
    p = [∇², linearOperator, mat1, mat2, r, a, divalphagrad]

    # Start progress logger if loggerFlag argument is 1
    loggerFlag==1 ? global_logger(TerminalLogger()) : nothing

    # Define split ODE problem
    prob = SplitODEProblem(DiffEqArrayOperator(linearOperator),splitNonlinearPart!,u0,(0.0,tMax),p)
    sol = solve(prob, ETDRK2(krylov=true, m=50), dt=δt, saveat=(tMax/100), rel_tol=0.001, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")

    # Calculate and plot free energy
    freeEnergies = freeEnergy(sol, ∇², mat1, mat2, nGrid, lSpace, r)

    if outputFlag==1
        params = @strdict nGrid lSpace r ϕ0 a δt tMax
        # Save variables and results to file
        fileName = savename(params, "jld2")
        @info "Saving data to output/$fileName"
        safesave("output/$fileName",@strdict sol freeEnergies params)
        # Plot results as animated gif
        if visualiseFlag==1
            visualise(sol, freeEnergies, params,"output/$fileName")
        end
    end

    return 1

end

export phaseFieldCrystal

end
