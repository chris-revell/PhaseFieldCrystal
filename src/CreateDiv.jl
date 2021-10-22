#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# 

module CreateDiv

# Import Julia packages
using LinearAlgebra
using SparseArrays
using LoopVectorization

arrayLoop(a,nGrid) = (nGrid+a-1)%(nGrid)+1

function createDiv(nGrid,h)

    adjX = spzeros(nGrid^2,nGrid^2)
    for x=1:nGrid
        for y=1:nGrid
            index1 = (x-1)*nGrid+y
            xNew = arrayLoop(x-1,nGrid)
            yNew = y
            index2 = (xNew-1)*nGrid + yNew
            adjX[index1,index2] = -1        
        end
    end

    adjY = spzeros(nGrid^2,nGrid^2)        
    for x=1:nGrid
        for y=1:nGrid
            index1 = (x-1)*nGrid+y # Index of grid point (x,y) when 2D array is flattened to a 1D vector            
            xNew = x
            yNew = arrayLoop(y-1,nGrid)
            index2 = (xNew-1)*nGrid + yNew
            adjY[index1,index2] = -1
        end
    end

    divX = (I+adjX)./h
    divY = (I+adjY)./h

    return divX, divY

end

export createDiv

end
