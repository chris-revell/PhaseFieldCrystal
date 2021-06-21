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

@inline @views function simulate(L, N ,r ,ϕ₀ ,α₀ ,q ,tMax)

    # Input parameters
    # L     Spatial dimensions of domain             (= 200.0 )
    # N     Number of grid points in each dimension  (= 200   )
    # r     Parameter in Swift-Hohenberg equation    (= -0.9  )
    # ϕ₀    Parameter in Swift-Hohenberg equation    (= -0.516)
    # α     Diffusivity                              (= 1.0   )
    # q     Parameter in Swift-Hohenberg equation    (= 0.1   )
    # tMax  Run time of simulation                   (= 2000.0)

    # Derived parameters
    h      = L/(N-1)    # Spatial separation of grid points
    Δ      = r+3.0*ϕ₀^2 # By definition
    outInt = tMax/100.0 # Output interval
    tspan  = (0.0,tMax) # Time span for solution

    # Create output folder and data files
    folderName = createRunDirectory(N,h,r,ϕ₀,Δ,α₀,q,outInt,tMax)

    # Set initial conditions: define arrays for calculations and set initial u0 order parameter field
    u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ = initialConditions(N,L,α₀)

    # Array of parameters to pass to solver
    p = [deriv, part1, part2, N, h, αᵢ, αⱼ, Δ, ϕ₀, q, graduᵢ, graduⱼ]

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(PFC!, u0, tspan, p)

    # Solve problem
    sol = solve(prob, Tsit5(), reltol=1e-6, saveat=outInt, maxiters=1e9)

    # Identify maximum and minimum values for colormap
    uMax = maximum(maximum.(sol.u))
    uMin = minimum(minimum.(sol.u))

    # Plot results as animated gif
    anim = @animate for i=1:(size(sol.t)[1])
       heatmap(sol.u[i][4:N+3,4:N+3].+ϕ₀,clims=(uMin+ϕ₀,uMax+ϕ₀),aspect_ratio=:equal,border=:none)
    end every 5
    gif(anim,"output/$folderName/anim.gif",fps=5)

    # Integrate over domain to check for mass conservation
    x = range(0, 1, length=N)
    y = range(0, 1, length=N)
    for u in sol.u
        display(integrate((x,y), u[3:N+2,3:N+2]))
    end

    return 1

end

export simulate

end
