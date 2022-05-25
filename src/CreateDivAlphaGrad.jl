#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
#

# module CreateDivAlphaGrad

# Import Julia packages
using LinearAlgebra
using SparseArrays

function createDivAlphaGrad(nX, nY, h, α)

    incidence = spzeros(2*nX*nY,nX*nY)
    for x=1:nX
        for y=1:nY
            # Index of grid point (x,y) when 2D array is flattened to a 1D vector
            indexVertex = (x-1)*nY+y

            # Loop over all edges neighbouring vertex (x,y)
            # Edges stored as vector of size 2*nX*nY.
            # First half of vector stores x directed edges; second half stores y directed edges.
            # Incidence matrix maps each vertex to 2 x directed edges with + and - directions, and 2 y directed edges with + and - directions

            # x dimension edge neighbours
            xNew = arrayLoop(x-1,nX)
            yNew = y
            indexEdge = (xNew-1)*nY + yNew
            incidence[indexEdge,indexVertex] = -1

            xNew = x
            yNew = y
            indexEdge = (xNew-1)*nY + yNew
            incidence[indexEdge,indexVertex] = 1

            # y dimension edge neighbours
            xNew = x
            yNew = arrayLoop(y-1,nY)
            indexEdge = (xNew-1)*nY + yNew + nX*nY
            incidence[indexEdge,indexVertex] = -1

            xNew = x
            yNew = y
            indexEdge = (xNew-1)*nY + yNew + nX*nY
            incidence[indexEdge,indexVertex] = 1

        end
    end
    incidence .= incidence./h

    incidenceTranspose = sparse(transpose(incidence))

    operator = -incidenceTranspose*α*incidence

    return operator

end

# export createDivAlphaGrad

# end
