#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
#

module CreateDivAlphaGrad

# Import Julia packages
using LinearAlgebra
using SparseArrays
using FromFile: @from

@from "ArrayLoop.jl" using ArrayLoop

function createDivAlphaGrad(nX, nY, h, α)

    incidence = spzeros(2*nX*nY,nX*nY) # Incidence matric mapping edges between grid points to grid points themselves 
    for x=1:nX
        for y=1:nY
            # Index of grid point (x,y) when 2D array is flattened to a 1D vector
            indexVertex = (x-1)*nY+y

            # Loop over all edges neighbouring vertex (x,y)
            # Edges stored as vector of size 2*nX*nY.
            # First half of vector stores x directed edges; second half stores y directed edges.
            # Incidence matrix maps each vertex to 2 x-directed edges with + and - directions, and 2 y-directed edges with + and - directions

            #  Find indices of x-direction edges surrounding grid point (x,y) ie (x-1,y) and (x,y)
            # Trailing edge:
            xNew = arrayLoop(x-1,nX)
            yNew = y
            indexEdge = (xNew-1)*nY + yNew        #  Convert cartesian indices of edges to within flattened vector
            incidence[indexEdge,indexVertex] = -1 #  Set incidence matrix component to -1 for trailing edge
            # Leading edge:
            xNew = x
            yNew = y
            indexEdge = (xNew-1)*nY + yNew        #  Convert cartesian indices of edges to within flattened vector
            incidence[indexEdge,indexVertex] = 1  #  Set incidence matrix component to 1 for leading edge

            #  Find indices of y-direction edges surrounding grid point (x,y) ie (x,y-1) and (x,y)
            #  Trailing edge
            xNew = x
            yNew = arrayLoop(y-1,nY)
            indexEdge = (xNew-1)*nY + yNew + nX*nY #  Convert cartesian indices of edges to within flattened vector, NB adding NxNy to put y-directed components in 2nd half of vector 
            incidence[indexEdge,indexVertex] = -1  #  Set incidence matrix component to -1 for trailing edge
            #  Leading edge
            xNew = x
            yNew = y
            indexEdge = (xNew-1)*nY + yNew + nX*nY #  Convert cartesian indices of edges to within flattened vector, NB adding NxNy to put y-directed components in 2nd half of vector 
            incidence[indexEdge,indexVertex] = 1   #  Set incidence matrix component to 1 for trailing edge

        end
    end
    incidence .= incidence./h   

    incidenceTranspose = sparse(transpose(incidence))

    operator = -incidenceTranspose*α*incidence

    return operator

end

export createDivAlphaGrad

end
