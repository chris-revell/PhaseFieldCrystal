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

@inline function boundaryConditions!(u, N, h)

    # Set 1st derivative at boundary to be zero
    u[:,3] .= u[:,5]
    u[:,N+4] .= u[:,N+2]
    u[3,:] .= u[5,:]
    u[N+4,:] .= u[N+2,:]

    # # Set 3rd derivative at boundary to be zero
    u[:,2]   .= u[:,6]
    u[2,:]   .= u[6,:]
    u[:,N+5] .= u[:,N+1]
    u[N+5,:] .= u[N+1,:]

    u[:,1]   .= u[:,7]
    u[1,:]   .= u[7,:]
    u[:,N+6] .= u[:,N]
    u[N+6,:] .= u[N,:]

    return nothing

end

export boundaryConditions!

end
