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

# ϕ̃ₜ = α∇²(Δϕ̃ + (q+∇²)²ϕ̃ + 3ϕ₀ϕ̃² + ϕ̃³)
# ϕ̃ₜ = α∇²(Δϕ̃ + ∇²(∇²ϕ̃ + qϕ̃) + q² + 3ϕ₀ϕ̃² + ϕ̃³)
# Δ = r+3ϕ₀²
# ϕ̃ = ϕ - ϕ₀
# From Archer, Knobloch 2012 Appendix B

@inline @views function PFC!(du, u, p, t)

    deriv, part1, part2, N, h, αᵢ, αⱼ, Δ, ϕ₀, q, graduᵢ, graduⱼ = p

    boundaryConditions!(u,N,h)

    ∇²!(deriv,u,N,h,0)

    part1 .= deriv .+ q.*u

    ∇²!(deriv, part1, N, h, 1)

    part2 .= deriv .+ Δ.*u .+ u.^3 .+ q^2 # .+ 3.0*ϕ₀.*u.^2

    grad!(graduᵢ, graduⱼ, part2, N, h)
    graduᵢ .*= αᵢ
    graduⱼ .*= αⱼ
    div!(du, graduᵢ, graduⱼ, N, h)

    return nothing

end

export PFC!

end
