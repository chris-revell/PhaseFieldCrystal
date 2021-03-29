#
#  Laplacian.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 29/03/2021.
#
#
# In place Laplacian function

module Laplacian

using LinearAlgebra

@inline @views function ∇²!(∇²u, u, N, h)

    for i = 2:N+1
        for j = 2:N+1
            ∇²u[i,j] = (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])/h^2
        end
    end

    return nothing

end

export ∇²!

end
