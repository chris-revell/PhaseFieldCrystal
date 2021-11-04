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

function createGrad(nGrid, h, α)

    incidence = spzeros(2*nGrid^2,nGrid^2)
    dx = [ 0, 0, 1, -1]
    dy = [-1, 1, 0,  0]
    for x=1:nGrid
        for y=1:nGrid
            indexVertex = (x-1)*nGrid+y # Index of grid point (x,y) when 2D array is flattened to a 1D vector

            # Loop over all edges neighbouring vertex (x,y)

            # x dimension edge neighbours
            xNew = arrayLoop(x-1,nGrid)
            yNew = arrayLoop(y,nGrid)
            indexEdge = (xNew-1)*nGrid + yNew
            incidence[indexEdge,indexVertex] = -1

            xNew = arrayLoop(x,nGrid)
            yNew = arrayLoop(y,nGrid)
            indexEdge = (xNew-1)*nGrid + yNew
            incidence[indexEdge,indexVertex] = 1

            # y dimension edge neighbours
            xNew = arrayLoop(x,nGrid)
            yNew = arrayLoop(y-1,nGrid)
            indexEdge = (xNew-1)*nGrid + yNew + nGrid^2
            incidence[indexEdge,indexVertex] = -1

            xNew = arrayLoop(x,nGrid)
            yNew = arrayLoop(y,nGrid)
            indexEdge = (xNew-1)*nGrid + yNew + nGrid^2
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
