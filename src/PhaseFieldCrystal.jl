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
using Model
using CreateRunDirectory
using InitialConditions
using Visualise
using ImportImage
using FreeEnergy

@inline @views function phaseFieldCrystal(imagePath,L,r,ϕ₀,α₀,q,tMax,visualiseFlag)

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
    h      = L/N        # Spatial separation of grid points
    outInt = tMax/100   # Output interval
    q²     = q^2        # Precalculate powers of q to reduce later calculations
    q⁴     = q^4        # Precalculate powers of q to reduce later calculations

    # Create output folder and data files
    folderName = createRunDirectory(L,N,h,r,ϕ₀,α₀,q,outInt,tMax)

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ = initialConditions(imageMask,L,N,α₀,ϕ₀,q)

    # Array of parameters to pass to solver
    p = [deriv, part1, part2, N, h, αᵢ, αⱼ, r, q, q², q⁴, graduᵢ, graduⱼ]

    ## Define ODE problem using discretised derivatives
    #prob = ODEProblem(PFC!, u0, tspan, p)

    # Start progress logger
    #global_logger(TerminalLogger())


    # Physics
    ttot   = 1.0          # total simulation time
    dt     = 0.2          # physical time step
    # Numerics
    nx, ny = N, N     # number of grid points
    tol    = 1e-6         # tolerance
    itMax  = 1e5          # max number of iterations
    # Derived numerics
    dx = h
    dy = h # grid size

    # Array allocation
    # (=graduⱼ) qHx    = zeros(nx-1, ny-2) # on staggered grid
    # (=graduᵢ) qHy    = zeros(nx-2, ny-1) # on staggered grid
    ResH   = zeros(nx, ny) # normal grid, without boundary points
    dudτ   = zeros(nx, ny) # normal grid, without boundary points
    dτ     = zeros(nx, ny) # normal grid, without boundary points
    uOld   = copy(u)
    dudt   = zeros(nx+6, ny+6)
    t      = 0.0
    it     = 0
    ittot  = 0

    # Physical time loop
    while t<ttot
        iter = 0
        err = 2*tol
        # Picard-type iteration
        while err>tol && iter<itMax

            dudt .= PFC!(dudt, u, p, t)

            ResH   .= inn(dudt) .- (inn(u) .- inn(uOld))/dt # residual of the PDE

            dHdtau .= ResH                                             # rate of change

            dtau   .= (1.0./(min(dx, dy)^2 ./inn(u).^npow./4.1) .+ 1.0/dt).^-1  # time step (obeys ~CFL condition)

            inn(u) .= inn(u) .+ dτ.*dHdτ              # update rule, sets the BC as u[1]=u[end]=0

            iter += 1

            err = norm(ResH)/length(ResH)
        end
        ittot += iter; it += 1; t += dt
        uOld .= u
    end



    # Calculate and plot free energy
    # freeEnergies = freeEnergy(sol, N, L, q, r, h)
    #
    # # Save variables and results to file
    # @info "Saving data to output/$folderName/data.jld2"
    # jldsave("output/$folderName/data.jld2";sol,imageMask,freeEnergies,N,L,q,r,h,folderName)
    #
    # # Plot results as animated gif
    # if visualiseFlag==1
    #     visualise(sol,N,h,folderName,freeEnergies,imageMask)
    # end

    return 1

end

@views   inn(A) = A[4:end-3,4:end-3]

export phaseFieldCrystal

end
