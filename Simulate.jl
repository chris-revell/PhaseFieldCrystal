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
using Plots
using Random
using NumericalIntegration

# Import local modules
using Model
using CreateRunDirectory
using InitialConditions

@inline @views function simulate(L,N,r,ϕ₀,α₀,q,tMax)

    # Input parameters
    # L     Spatial dimensions of domain             (= 200.0 )
    # N     Number of grid points in each dimension  (= 200   )
    # r     Parameter in Swift-Hohenberg equation    (= -0.9  )
    # ϕ₀    Mean order parameter across domain       (= -0.516)
    # α     Diffusivity                              (= 1.0   )
    # q     Parameter in Swift-Hohenberg equation    (= 0.1   )
    # nPlot Number of plots produced by return       (= 100   )
    # tMax  Run time of simulation                   (= 2000.0)

    # Derived parameters
    h      = L/(N-1)    # Spatial separation of grid points
    outInt = tMax/100 # Output interval
    tspan  = (0.0,tMax) # Time span for solution

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ,ϕ₀Real = initialConditions(L,N,α₀,ϕ₀)

    # Create output folder and data files
    folderName = createRunDirectory(L,N,h,r,ϕ₀,α₀,q,outInt,tMax,ϕ₀Real)

    # Array of parameters to pass to solver
    p = [deriv, part1, part2, N, h, αᵢ, αⱼ, r, q, graduᵢ, graduⱼ]

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(PFC!, u0, tspan, p)

    # Solve problem
    sol = solve(prob, Tsit5(), reltol=1e-2, saveat=outInt, maxiters=1e9)

    # # Identify maximum and minimum values for colormap
    # uMax = maximum(maximum.(sol.u[2:end]))
    # uMin = minimum(minimum.(sol.u[2:end]))

    # Plot results as animated gif
    anim = @animate for i=1:(size(sol.t)[1])
       heatmap(sol.u[i][4:N+3,4:N+3],clims=(-1,1),aspect_ratio=:equal,border=:none)
    end every 5
    gif(anim,"output/$folderName/anim.gif",fps=2)

    # Integrate over domain to check for mass conservation
    # x = range(0, L, length=N)
    # y = range(0, L, length=N)
    # for u in sol.u
    #     display(integrate((x,y), u[4:N+3,4:N+3])/(L*L))
    # end

    return 1

end

export simulate

end
