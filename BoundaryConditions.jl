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
    # Set 1st ghost points [2,:] and [:,2] to equal penultimate internal points.
    u[:,2] .= u[:,4]
    u[:,N+3] .= u[:,N+1]
    u[2,:] .= u[4,:]
    u[N+3,:] .= u[N+1,:]

    # # Set 3rd derivative at boundary to be zero
     u[:,1]   .= u[:,5]
     u[1,:]   .= u[5,:]
     u[:,N+4] .= u[:,N]
     u[N+4,:] .= u[N,:]


     # u[2,1]   = u[2,2]
     # u[1,2]   = u[2,2]
     # u[N+3,1] = u[N+3,2]
     # u[N+4,2] = u[N+3,2]
     # u[1,N+3] = u[2,N+3]
     # u[2,N+4] = u[2,N+3]
     # u[N+4,N+3] = u[N+3,N+3]
     # u[N+3,N+4] = u[N+3,N+3]



    return nothing

end

export boundaryConditions!

end
