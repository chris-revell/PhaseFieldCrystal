#
#  BoundaryConditions.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# In place operations to ensure no flux boundary condition in heat equation

module BoundaryConditions

# Import Julia packages
using LinearAlgebra

# Import local modules
# include("<Module>.jl"); using .Module

@views function boundaryConditions!(u, nGrid, inset)

    uTmp = reshape(u,(nGrid+6,nGrid+6))

    # Set 1st derivative at boundary to be zero by reflecting values around edge point
    uTmp[:,inset]   .= uTmp[:,inset+2]
    uTmp[:,nGrid+inset+1] .= uTmp[:,nGrid+inset-1]
    uTmp[inset,:]   .= uTmp[inset+2,:]
    uTmp[nGrid+inset+1,:] .= uTmp[nGrid+inset-1,:]

    # Set 3rd derivative at boundary to be zero by reflecting values around edge point
    uTmp[:,inset-1]   .= uTmp[:,inset+3]
    uTmp[inset-1,:]   .= uTmp[inset+3,:]
    uTmp[:,nGrid+inset+2] .= uTmp[:,nGrid+inset-2]
    uTmp[nGrid+inset+2,:] .= uTmp[nGrid+inset-2,:]

    # Set 5th derivative at boundary to be zero by reflecting values around edge point
    uTmp[:,inset-2]   .= uTmp[:,inset+4]
    uTmp[inset-2,:]   .= uTmp[inset+4,:]
    uTmp[:,nGrid+inset+3] .= uTmp[:,nGrid+inset-3]
    uTmp[nGrid+inset+3,:] .= uTmp[nGrid+inset-3,:]    

    return nothing

end

export boundaryConditions!

end
