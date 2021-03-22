#
#  Simulate.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#

module Model

using SparseArrays

# Function for time derivative of Cahn-Hilliard using discretised derivatives
@inline @views function dϕ!(du, u, p, t)

    Δ₁, Q, α, γ = p # Opening up parameters array argument to obtain: 2nd order derivative operator, boundary condition operator, diffusion coefficient, gamma coefficient

    # Two components of C-H equation
    ∇²u = Δ₁*Q*u .+ (Δ₁*Q*u')'
    part1 = u.^3 .- u .- γ.*∇²u
    # Combining the above components to obtain du/dt
    du .= α.*(Δ₁*Q*part1 .+ (Δ₁*Q*part1')')

end

export dϕ!

end
