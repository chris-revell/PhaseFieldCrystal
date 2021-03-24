#
#  Simulate.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#

module Simulate

# Import Julia packages
using DiffEqOperators
using DifferentialEquations
using LinearAlgebra
using Plots
#using JLD2

# Import local modules
using Model
using CreateRunDirectory

@inline @views function simulate()

    L      = 20      # Number of grid pointsin each dimension
    dx     = 1/(L-1) # Spatial separation of grid points
    r      = 1.0     # Δ in paper renamed to r to avoid confusion with derivatives
    q      = 1.0
    ϕ₀     = 1.0
    tMax   = 0.0001
    outInt = tMax/100.0
    tspan  = (0.0,tMax)   # Time span for solution

    # Create output files
    folderName = createRunDirectory(L,dx,r,q,ϕ₀,tMax,outInt)

    # Random initial distribution
    u0             = rand(L,L)
    uᵀ             = zeros(L,L)
    firstDimTerm   = zeros(L,L)
    secondDimTerm  = zeros(L,L)
    secondDimTermᵀ = zeros(L,L)

    # Dirichlet boundary condition
    Q = Neumann0BC(dx,3)

    # Spatial derivative operators
    ∇₂ = CenteredDifference{1}(2, 3, dx, L)
    ∇₄ = CenteredDifference{1}(4, 3, dx, L)
    ∇₆ = CenteredDifference{1}(6, 3, dx, L)

    p = [∇₂,∇₄,∇₆,Q,r,q,ϕ₀,uᵀ,firstDimTerm,secondDimTerm,secondDimTermᵀ]      # Array of parameters to pass to solver

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(dϕ!, u0, tspan, p)

    # Solve problem
    sol = solve(prob, Tsit5(), reltol=1e-6, saveat=outInt, maxiters=1e7)

    # Plot results as animated gif
    anim = @animate for i=1:(size(sol.t)[1])
       heatmap(sol.u[i],clims=(-1,1),aspect_ratio=:equal)
       #surface(sol.u[i],zlims=(0,1))
    end every 5
    gif(anim,"output/$folderName/anim.gif",fps=2)

    return 1

end

export simulate

end
