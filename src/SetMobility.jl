#
#  SetMobility.xl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 08/11/2021.
#
#
# Function to set mobility matrix based on imported image


# module SetMobility

using LinearAlgebra
using SparseArrays

function setMobility(nX,nY,imageMask)

    # Mobility stored as a vector twice as long as system vector
    # 2 mobility values for each grid point corresponding to edges in x and y directions.
    αVec = zeros(2*nX*nY)
    αᵢTmp = ones(nY,nX)
    αⱼTmp = ones(nY,nX)
    for x=1:nX
        for y=1:nY
            if imageMask[y,x] == 0.0 || imageMask[(nY+y+1-1)%(nY)+1,x] == 0.0
                αᵢTmp[y,x] = 0.0
            end
        end
    end
    for x=1:nX
        for y=1:nY
            if imageMask[y,x] == 0.0 || imageMask[y,(nX+x+1-1)%(nX)+1] == 0.0
                αⱼTmp[y,x] = 0.0
            end
        end
    end
    αVec[1:nX*nY] .= reshape(αᵢTmp,nX*nY)
    αVec[1+nX*nY:end] .= reshape(αⱼTmp,nX*nY)
    α = spdiagm(αVec)

    return α

end

# export setMobility

# end
