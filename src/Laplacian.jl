#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# In place Laplacian function

module Laplacian

# Import Julia packages
using LinearAlgebra
using LoopVectorization
using SparseArrays

# Import local modules
# include("<Module>.jl"); using .Module

function createLaplacian(N, h)

    adj = spzeros(N^2,N^2)
    dx = [ 0, 0, 1, -1]
    dy = [-1, 1, 0,  0]
    for x=1:N
        for y=1:N
            index = (x-1)*N+y
            for ne=1:length(dx)
                newx = x + dx[ne]
                newy = y + dy[ne]
                if newx>0 && newx <=N && newy>0 && newy<=N
                    index2 = (newx-1)*N + newy
                    adj[index,index2] = 1
                end
            end
        end
    end

    degree = spdiagm(0=>sum(adj, dims=2)[:,1])
    ∇² = (degree-adj)./h^2

    return ∇²

end

export createLaplacian

end
