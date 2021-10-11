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

# Import local modules
# include("<Module>.jl"); using .Module

arrayLoop(a,nGrid) = (nGrid+a-1)%(nGrid)+1

function createDiv

    adj = spzeros(nGrid^2,nGrid^2)
    # List of neighbour indices in 
    dx = [ 0, 0, 1, -1]
    dy = [-1, 1, 0,  0]
    for x=1:nGrid
        for y=1:nGrid
            index1 = (x-1)*nGrid+y # Index of grid point (x,y) when 2D array is flattened to a 1D vector 
            # Loop over all neighbours of (x,y)
            for i=1:length(dx)
                xNew = arrayLoop(x+dx[i],nGrid) # Find (x,y) indices of neighbouring grid point, introducing periodicity with arrayLoop
                yNew = arrayLoop(y+dy[i],nGrid) # Find (x,y) indices of neighbouring grid point, introducing periodicity with arrayLoop
                index2 = (xNew-1)*nGrid + yNew
                adj[index1,index2] = 1                
            end
        end
    end

    degree = spdiagm(0=>sum(adj, dims=2)[:,1])
    ∇² = (adj-degree)./h^2

    return ∇²

end

export createDiv

end
