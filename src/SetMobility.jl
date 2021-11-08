#
#  SetMobility.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 08/11/2021.
#
#
# Function to set mobility matrix based on imported image


module SetMobility

using LinearAlgebra
using SparseArrays

function setMobility(nGrid,imageMask)

    αVec = zeros(2*nGrid^2)
    αᵢTmp = ones(nGrid,nGrid)
    αⱼTmp = ones(nGrid,nGrid)
    for j=1:nGrid
        for i=1:nGrid
            if imageMask[i,j] == 0.0 || imageMask[(nGrid+i+1-1)%(nGrid)+1,j] == 0.0
                αᵢTmp[i,j] = 0.0
            end
        end
    end
    for j=1:nGrid
        for i=1:nGrid
            if imageMask[i,j] == 0.0 || imageMask[i,(nGrid+j+1-1)%(nGrid)+1] == 0.0
                αⱼTmp[i,j] = 0.0
            end
        end
    end
    αVec[1:nGrid^2] .= reshape(αᵢTmp,nGrid^2)
    αVec[1+nGrid^2:end] .= reshape(αⱼTmp,nGrid^2)
    α = spdiagm(αVec)

    return α

end

export setMobility

end
