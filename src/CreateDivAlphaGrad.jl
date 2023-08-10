#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# Due to the spatially varying diffusivity, we must consider the outermost differential operators ($\nabla\cdot\alpha\nabla$) separately. In this instance, we construct a directed incidence matrix $\mathcal{A}$. Incidence matrix $\mathcal{A}$ is a $2NM\times NM$ matrix mapping edges between grid points (vertices) in rows to grid points in columns. Note that there are twice as many edges as vertices, hence why $\mathcal{I}$ has twice as many rows as columns.
# Incidence matrix $\mathcal{A}$ has a value of 1 for the components mapping grid point $(i,j)$ to leading edge $(i,j)$, which itself links grid points $(i,j)$ and $(i+1,j)$; a value of -1 for the component mapping grid point $(i,j)$ to trailing edge $(i-1,j)$, which links grid point $(i,j)$ and $(i-1,j)$, with a similar protocol for y-directed edges. With $\alpha$ a diagonal $2NM\times 2NM$ matrix containing diffusiivity values, the external Laplacian operator $\nabla\cdot\alpha\nabla = -\mathcal{I}^T\alpha\mathcal{I}$.


module CreateDivAlphaGrad

# Import Julia packages
using LinearAlgebra
using SparseArrays
using DrWatson
using FromFile: @from

@from "$(srcdir("ArrayLoop.jl"))" using ArrayLoop

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
