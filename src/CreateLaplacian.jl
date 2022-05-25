#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# In place Laplacian function

# module CreateLaplacian

# Import Julia packages
using LinearAlgebra
using SparseArrays

function createLaplacian(nX, nY, h)

    adj = spzeros(nX*nY,nX*nY)
    # List of neighbour indices in
    dx = [ 0, 0, 1, -1]
    dy = [-1, 1, 0,  0]
    for x=1:nX
        for y=1:nY
            index1 = (x-1)*nY+y # Index of grid point (x,y) when 2D array is flattened to a 1D vector
            # Loop over all neighbours of (x,y)
            for i=1:length(dx)
                xNew = arrayLoop(x+dx[i],nX) # Find (x,y) indices of neighbouring grid point, introducing periodicity with arrayLoop
                yNew = arrayLoop(y+dy[i],nY) # Find (x,y) indices of neighbouring grid point, introducing periodicity with arrayLoop
                index2 = (xNew-1)*nY + yNew
                adj[index1,index2] = 1
            end
        end
    end

    degree = spdiagm(0=>sum(adj, dims=2)[:,1])
    ∇² = (adj-degree)./h^2

    return ∇²

end

# export createLaplacian

# end
