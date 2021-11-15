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

function setMobility(nY,nX,imageMask)

    αVec = zeros(2*nY*nX)
    αᵢTmp = ones(nX,nY)
    αⱼTmp = ones(nX,nY)
    for j=1:nY
        for i=1:nX
            if imageMask[i,j] == 0.0 || imageMask[(nX+i+1-1)%(nX)+1,j] == 0.0
                αᵢTmp[i,j] = 0.0
            end
        end
    end
    for j=1:nY
        for i=1:nX
            if imageMask[i,j] == 0.0 || imageMask[i,(nY+j+1-1)%(nY)+1] == 0.0
                αⱼTmp[i,j] = 0.0
            end
        end
    end
    αVec[1:nY*nX] .= reshape(αᵢTmp,nY*nX)
    αVec[1+nY*nX:end] .= reshape(αⱼTmp,nY*nX)
    α = spdiagm(αVec)

    return α

end

export setMobility

end
