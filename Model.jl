#
#  Model.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#
# Function for time derivative of Cahn-Hilliard equation using discretised derivatives

module Model

# Julia packages
using LinearAlgebra

# Local modules
using BoundaryConditions
using Laplacian

@inline @views function CH!(du, u, p, t)

    ∇²u, N, h, α, γ, part1 = p

    boundaryConditions!(u,N,h)

    ∇²!(∇²u,u,N,h,0)

    part1 .= α.*(u.^3 .- u .- γ.*∇²u)

    ∇²!(du,part1,N,h,1)

end

export CH!

end
