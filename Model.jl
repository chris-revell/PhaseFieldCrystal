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

    ∇₂,∇₄,∇₆,Q,r,q,ϕ₀ = p # Opening up parameters array argument

    # ∇²u = ∇₂*Q*u .+ (∇₂*Q*u')'
    # part1 = u.^3 .- u .- ∇²u
    # du .= (∇₂*Q*part1 .+ (∇₂*Q*part1')')
    du .= (r.*∇₂*Q*u .+ q².*∇₂*Q*u .+ q².*∇₄*Q*u .+ ∇₆*Q*u .+ 3.0*ϕ₀.*∇₂*Q*(u.^2) .+ ∇₂*Q*(ϕ.^3)) .+ (r.*∇₂*Q*u' .+ q².*∇₂*Q*u' .+ q².*∇₄*Q*u' .+ ∇₆*Q*u' .+ 3.0*ϕ₀.*∇₂*Q*(u.^2)' .+ ∇₂*Q*(ϕ.^3)')'

end

export dϕ!

end
