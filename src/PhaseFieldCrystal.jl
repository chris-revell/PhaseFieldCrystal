#
#  PhaseFieldCrystal.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#
# Input parameters:
# imagePath      Location in filesystem of image to use as solution domain
# lX             Spatial dimensions of domain
# r              Parameter in Swift-Hohenberg equation     (eg. = 0.5  )
# ϕ0             Mean order parameter across domain        (eg. = -0.37)
# a              Parameter in splitting scheme             (eg. = 2.0  )
# λ              Lengthscale of Gaussian Random Field      (eg. = 1.0  )
# δt             Time step for implicit semi-linear scheme (eg. = 0.01 )
# tMax           Run time of simulation                    (eg. = 20.0 )
# outCount       Number of data output timepoints          (eg. = 100  )
# loggerFlag     Display progress bar in REPL when =1; Always set to 0 on HPC
# outputFlag     Flag to control whether data are saved to file (=1 or 0)
# visualiseFlag  Flag to control whether data are plotted  (=1 or 0)
# freeEnergyFlag Flag to control whether free energues are calculated  (=1 or 0)
# nBlasThreads   Number of threads to allow for BLAS operations
# subFolder      Optional parameter to specify subfolder within data/sims in which to store results 

# Derived parameters:
# imageMask     Binarised array of 0s and 1s representing inter-cell spaces and cell regions in image
# nY            Number of grid points in x dimension (horizontal in image; 2nd matrix dimension)
# nX            Number of grid points in y dimension (vertical in image; 1st matrix dimension)
# h             Spatial step between grid points

module PhaseFieldCrystal

# Import Julia packages
using DifferentialEquations
using LinearAlgebra
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using SparseArrays
using DrWatson
using Dates
using Base.Filesystem
using FromFile: @from

# Import local files
@from "Model.jl" using Model
@from "CreateLaplacian.jl" using CreateLaplacian
@from "CreateDivAlphaGrad.jl" using CreateDivAlphaGrad
@from "InitialConditions.jl" using InitialConditions
@from "Visualise.jl" using Visualise
@from "FreeEnergy.jl" using FreeEnergy
@from "ImportImage.jl" using ImportImage
@from "SetMobility.jl" using SetMobility

function phaseFieldCrystal(imagePath,lX,r,ϕ0,m,a,λ,δt,tMax,outCount,loggerFlag,outputFlag,visualiseFlag,freeEnergyFlag,nBlasThreads;subFolder="")

    BLAS.set_num_threads(nBlasThreads)

    imageMask, nX, nY = importImage(imagePath)

    α = setMobility(nX,nY,imageMask)

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,mat1,mat2,h = initialConditions(imageMask,lX,nX,nY,ϕ0,λ,m)

    # Create finite difference matrices for given system parameters
    ∇² = createLaplacian(nX,nY,h)
    divalphagrad = createDivAlphaGrad(nX,nY,h,α)

    # Create matrix for linear component of PFC equation
    linearOperator = divalphagrad.*(1.0-r+a) .+ divalphagrad*∇²*∇²

    # Array of parameters to pass to solver
    p = [∇², linearOperator, mat1, mat2, r, a, divalphagrad]

    # Start progress logger if loggerFlag argument is 1
    loggerFlag==1 ? global_logger(TerminalLogger()) : nothing

    # Define split ODE problem
    prob = SplitODEProblem(DiffEqArrayOperator(linearOperator),splitNonlinearPart!,u0,(0.0,tMax),p)     #,tstops=tStopsArray)
    sol = solve(prob, ETDRK2(krylov=true, m=50), dt=δt, saveat=tMax/outCount, reltol=0.001, progress=(loggerFlag==1), progress_steps=10, progress_name="PFC model")

    # Calculate free energy at each time point of solution
    freeEnergyFlag==1 ? freeEnergies = freeEnergy(sol, ∇², mat1, mat2, nX, nY, lX, r) : nothing

    # Save results in JLD2 format with unique filename
    if outputFlag==1
        maskFileName = splitpath(imagePath)[end][1:end-4]
        mkpath(datadir("sims",subFolder))
        params = @strdict  maskFileName ϕ0 r m λ nX nY lX a δt tMax
        # Create filename from parameters; prefix filename with current data and time
        fileName = savename(params,connector="",ignores=["a"])
        # Save variables and results to file
        u = sol.u        
        t = sol.t
        safesave(datadir("sims",subFolder,"$fileName.jld2"),@strdict u t ϕ0 r m λ nX nY lX a δt tMax)
        # Plot results as animated gif and free energies as png
        if (visualiseFlag==1)
            visualise(u, t, ϕ0, r, m, nX, nY, lX, a, δt, tMax, subFolder, fileName, freeEnergyFlag)
        end
        @info "Saved data to $(datadir("sims",subFolder,"$fileName.jld2"))"
    end

    return nothing

end

export phaseFieldCrystal

end

# Block for adjusting timestep during run
    # tStopsArray = Float64[]
    # push!(tStopsArray,0.0)
    # tTmp = 0.0
    # while tTmp<100.0
    #     tTmp += δt
    #     push!(tStopsArray,tTmp)
    # end
    # while tTmp<tMax
    #     tTmp += 2.0*δt
    #     push!(tStopsArray,tTmp)
    # end