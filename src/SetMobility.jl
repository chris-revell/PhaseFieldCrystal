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

function setMobility(nX,nY,imageMask)

    # Mobility stored as a vector twice as long as system vector
    # 2 mobility values for each grid point corresponding to edges in x and y directions.
    αVec = zeros(2*nX*nY)
    αᵢTmp = ones(nY,nX)
    αⱼTmp = ones(nY,nX)
    for j=1:nX
        for i=1:nY
            if imageMask[i,j] == 0.0 || imageMask[(nY+i+1-1)%(nY)+1,j] == 0.0
                αᵢTmp[i,j] = 0.0
            end
        end
    end
    for j=1:nX
        for i=1:nY
            if imageMask[i,j] == 0.0 || imageMask[i,(nX+j+1-1)%(nX)+1] == 0.0
                αⱼTmp[i,j] = 0.0
            end
        end
    end
    αVec[1:nX*nY] .= reshape(αᵢTmp,nX*nY)
    αVec[1+nX*nY:end] .= reshape(αⱼTmp,nX*nY)
    α = spdiagm(αVec)

    return α

end

export setMobility

end
