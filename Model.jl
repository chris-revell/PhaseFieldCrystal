#
#  Model.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#

module Model

# Julia packages
using LinearAlgebra

# Local modules
using BoundaryConditions
using Laplacian
using Grad
using Div

# ϕ̃ₜ = α∇²(rϕ + (q²+∇²)²ϕ + ϕ³)
# ϕ̃ₜ = α∇²(rϕ + ∇²(∇²ϕ + q²ϕ) + q⁴ + ϕ³)

# From Thiele, Knobloch 2013 Equation 3

@inline @views function PFC!(du, u, p, t)

    deriv, part1, part2, N, h, αᵢ, αⱼ, r, q, graduᵢ, graduⱼ = p

    boundaryConditions!(u,N,h)

    ∇²!(deriv,u,N,h,0)

    part1 .= deriv .+ (q^2).*u

    ∇²!(deriv, part1, N, h, 1)

    part2 .= r.*u .+ deriv .+ q^4 .+ u.^3

    grad!(graduᵢ, graduⱼ, part2, N, h)
    graduᵢ .*= αᵢ
    graduⱼ .*= αⱼ
    div!(du, graduᵢ, graduⱼ, N, h)

    return nothing

end

export PFC!

end
