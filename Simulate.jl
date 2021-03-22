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
using Plots
using LinearAlgebra
using JLD2

# Import local modules
using Model
using CreateRunDirectory

@inline @views function simulate()

    L = 50              # Number of grid pointsin each dimension
    dx = 1/(L-1)        # Spatial separation of grid points
    α = 0.5             # Diffusion coefficient
    γ = 0.00005         # sqrt(γ) gives the length of the transition regions between domains
    tMax = 0.5
    outputInterval = tMax/100.0
    tspan = (0.0,tMax)   # Time span for solution

    # Create output files
    folderName = createRunDirectory(L,dx,α,γ,tMax)

    # Random initial distribution
    u0 = 2.0.*rand(Float64,L,L).-1.0

    # Dirichlet boundary condition
    #Q = Dirichlet0BC(eltype(u0))
    Q = Neumann0BC(dx,3)

    # Second spatial derivative operator
    Δ₁ = CenteredDifference{1}(2, 3, dx, L)

    p = [Δ₁,Q,α,γ]      # Array of parameters to pass to solver

    # Define ODE problem using discretised derivatives
    prob = ODEProblem(dϕ!, u0, tspan, p)
    # Solve problem
    sol = solve(prob, reltol=1e-6, saveat=outputInterval, maxiters=1e7)

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
