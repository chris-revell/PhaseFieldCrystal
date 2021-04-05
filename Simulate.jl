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

@inline @views function simulate()

    L      = 50.0
    N      = 50         # Number of grid pointsin each dimension
    h      = L/(N-1)    # Spatial separation of grid points
    r      = -0.9
    ϕ₀     = -0.516
    Δ      = r+3.0*ϕ₀^2 # By definition
    α      = 1.0        #
    q      = 0.1        #
    tMax   = 100.0      #
    outInt = tMax/100.0 # Output interval
    tspan  = (0.0,tMax) # Time span for solution

    # Create output files
    folderName = createRunDirectory(N,h,r,ϕ₀,Δ,α,q,outInt,tMax)

    # Random initial distribution
    u0    = (2.0.*rand(N+6,N+6).-1.0)
    deriv = zeros(N+6,N+6)
    part1 = zeros(N+6,N+6)
    part2 = zeros(N+6,N+6)

    # Array of parameters to pass to solver
    p = [deriv, part1, part2, N, h, α, Δ, ϕ₀, q]

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(PFC!, u0, tspan, p)

    # Solve problem
    sol = solve(prob, Tsit5(), reltol=1e-6, saveat=outInt, maxiters=1e9)

    # Plot results as animated gif
    anim = @animate for i=1:(size(sol.t)[1])
       heatmap(sol.u[i][4:N+3,4:N+3],clims=(-1.0,1.0),aspect_ratio=:equal,legend=:false,border=:none)
    end every 5
    gif(anim,"output/$folderName/anim.gif",fps=5)

    # x = range(0, 1, length=N)
    # y = range(0, 1, length=N)
    # for u in sol.u
    #     display(integrate((x,y), u[3:N+2,3:N+2]))
    # end

    return 1

end

export simulate

end
