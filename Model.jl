#
#  Model.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#
# Function for time derivative of heat equation using discretised derivatives

module Model

# Julia packages
using LinearAlgebra
using Laplacian

# Local modules
using BoundaryConditions
using Laplacian

@inline @views function CH!(du, u, p, t)

    ∇²u, N, h, α = p

    boundaryConditions!(u,N)

    ∇²!(∇²u,u,N,h)

    du .= α.*∇²u

end

export CH!

end
