#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
#

module CreateGrad

# Import Julia packages
using LinearAlgebra
using SparseArrays
using LoopVectorization

arrayLoop(a,nGrid) = (nGrid+a-1)%(nGrid)+1

function createGrad(nY, nX, h, α)

    incidence = spzeros(2*nY*nX,nY*nX)
    dx = [ 0, 0, 1, -1]
    dy = [-1, 1, 0,  0]
    for x=1:nY
        for y=1:nX
            # Index of grid point (x,y) when 2D array is flattened to a 1D vector
            indexVertex = (x-1)*nX+y

            # Loop over all edges neighbouring vertex (x,y)

            # x dimension edge neighbours
            xNew = arrayLoop(x-1,nY)
            yNew = arrayLoop(y,nX)
            indexEdge = (xNew-1)*nX + yNew
            incidence[indexEdge,indexVertex] = -1

            xNew = arrayLoop(x,nY)
            yNew = arrayLoop(y,nX)
            indexEdge = (xNew-1)*nX + yNew
            incidence[indexEdge,indexVertex] = 1

            # y dimension edge neighbours
            xNew = arrayLoop(x,nY)
            yNew = arrayLoop(y-1,nX)
            indexEdge = (xNew-1)*nX + yNew + nY*nX
            incidence[indexEdge,indexVertex] = -1

            xNew = arrayLoop(x,nY)
            yNew = arrayLoop(y,nX)
            indexEdge = (xNew-1)*nX + yNew + nY*nX
            incidence[indexEdge,indexVertex] = 1

        end
    end
    incidence .= incidence./h

    incidenceTranspose = sparse(transpose(incidence))

    operator = -incidenceTranspose*α*incidence

    return operator

end

export createGrad

end
