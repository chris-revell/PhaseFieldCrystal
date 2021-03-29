#
#  BoundaryConditions.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# In place operations to ensure no flux boundary condition in heat equation

module BoundaryConditions

using LinearAlgebra

@inline @views function boundaryConditions!(u, N)

    u[:,1] .= u[:,3]
    u[:,N+2] .= u[:,N]
    u[1,:] .= u[3,:]
    u[N+2,:] .= u[N,:]

    return nothing

end

export boundaryConditions!

end
