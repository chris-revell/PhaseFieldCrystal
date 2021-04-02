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

@inline @views function simulate(N,α,γ,tMax)

    #N      = 100         # Number of grid pointsin each dimension
    h      = 1/(N-1)    # Spatial separation of grid points
    #α      = 0.5        # Diffusion coefficient
    #γ      = 0.0005    # sqrt(γ) gives the length of the transition regions between domains
    #tMax   = 5.0
    outInt = tMax/100.0 # Output interval
    tspan  = (0.0,tMax) # Time span for solution

    # Create output files
    folderName = createRunDirectory(N,h,α,γ,outInt,tMax)

    # Random initial distribution
    u0    = 2.0.*rand(N+4,N+4).-1.0
    ∇²u   = zeros(N+4,N+4)
    part1 = zeros(N+4,N+4)

    # Array of parameters to pass to solver
    p = [∇²u, N, h, α, γ, part1]

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(CH!, u0, tspan, p)

    # Solve problem
    sol = solve(prob, Tsit5(), reltol=1e-6, saveat=outInt, maxiters=1e7)

    # Plot results as animated gif
    anim = @animate for i=1:(size(sol.t)[1])
       heatmap(sol.u[i],clims=(-1.0,1.0),aspect_ratio=:equal,legend=:false,border=:none)
    end every 5
    gif(anim,"output/$folderName/anim.gif",fps=5)

    x = range(0, 1, length=N)
    y = range(0, 1, length=N)
    for u in sol.u
        display(integrate((x,y), u[3:N+2,3:N+2]))
    end

    return 1

end

export simulate

end
