#
#  Simulate.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#

module Model

using DiffEqOperators
using LinearAlgebra

# Function for conserved Swift-Hohenberg equation using discretised derivatives
@inline @views function dϕ!(du, u, p, t)
    # ϕ̇ = α∇²(rϕ + (q² + ∇²)²ϕ + 3ϕ₀ϕ² + ϕ³)
    # ϕ̇ = r∇²ϕ + ∇²q²ϕ + q²∇⁴ϕ + ∇⁶ϕ + 3ϕ₀∇²ϕ² + ∇²ϕ³

    ∇₂,∇₄,∇₆,Q,r,q,ϕ₀,uᵀ,firstDimTerm,secondDimTerm,secondDimTermᵀ = p # Opening up parameters array argument

    transpose!(uᵀ,u)

    firstDimTerm .= ∇₂*Q*(r.*u) .+ ∇₂*Q*(q^2.0 .*u) .+ ∇₄*Q*(q^2.0 .*u) .+ ∇₆*Q*u .+ ∇₂*Q*(3.0*ϕ₀.*u.^2) .+ ∇₂*Q*(u.^3)

    secondDimTerm .= ∇₂*Q*(r.*uᵀ) .+ ∇₂*Q*(q^2.0 .*uᵀ) .+ ∇₄*Q*(q^2.0 .*uᵀ) .+ ∇₆*Q*(uᵀ) .+ ∇₂*Q*(3.0*ϕ₀.*uᵀ.^2) .+ ∇₂*Q*(uᵀ.^3)

    transpose!(secondDimTermᵀ,secondDimTerm)

    du .=  firstDimTerm .+ secondDimTermᵀ

end

export dϕ!

end
