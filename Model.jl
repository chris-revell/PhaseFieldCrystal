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

    # Unpack parameter list
    deriv, part1, part2, N, h, αᵢ, αⱼ, r, q, q², q⁴, graduᵢ, graduⱼ = p

    # Set values of ghost points to ensure zero flux at boundary
    boundaryConditions!(u,N,h)

    # Find Laplacian of u
    ∇²!(deriv,u,N,h,0)

    # Calculate inner component (∇²ϕ + q²ϕ)
    part1 .= deriv .+ q².*u    

    # Find Laplacian of (∇²ϕ + q²ϕ)
    ∇²!(deriv, part1, N, h, 1)

    # Calculate full term within outermost Laplacian (rϕ + ∇²(∇²ϕ + q²ϕ) + q⁴ + ϕ³) = rϕ + ∇²(part1) + q⁴ + ϕ³
    part2 .= r.*u .+ deriv .+ q⁴ .+ u.^3    

    # Find grad of (rϕ + ∇²(∇²ϕ + q²ϕ) + q⁴ + ϕ³) and multiply by spatially varying diffusivity
    grad!(graduᵢ, graduⱼ, part2, αᵢ, αⱼ, N, h)

    # Find divergence of the result from the last term
    div!(du, graduᵢ, graduⱼ, N, h)

    return nothing

end

export PFC!

end
