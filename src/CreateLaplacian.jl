#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# In place Laplacian function

module CreateLaplacian

# Import Julia packages
using LinearAlgebra
using SparseArrays
using FromFile: @from

@from "ArrayLoop.jl" using ArrayLoop

function createLaplacian(nX, nY, h)

    adj = spzeros(nX*nY,nX*nY)
    # List of index steps needed to reach neighbours in 5-point Von Neumann stencil
    dx = [ 0, 0, 1, -1]
    dy = [-1, 1, 0,  0]
    for x=1:nX
        for y=1:nY
            flattenedIndex = (x-1)*nY+y # Index of grid point (x,y) when 2D array is flattened to a 1D vector
            # Loop over all neighbours of (x,y)
            for i=1:length(dx)
                xNeighbour = arrayLoop(x+dx[i],nX) # Find (x,y) indices of neighbouring grid point, introducing periodicity with arrayLoop
                yNeighbour = arrayLoop(y+dy[i],nY) # Find (x,y) indices of neighbouring grid point, introducing periodicity with arrayLoop
                flattenedIndexNeighbour = (xNeighbour-1)*nY + yNeighbour
                adj[flattenedIndex,flattenedIndexNeighbour] = 1 # Set corresponding component of adj to 1, indicating adjacency in 2D of grid points corresponding to flattenedIndex and flattenedIndexNeighbour
            end
        end
    end

    degree = spdiagm(0=>sum(adj, dims=2)[:,1])
    ∇² = (adj-degree)./h^2

    return ∇²

end

export createLaplacian

end
