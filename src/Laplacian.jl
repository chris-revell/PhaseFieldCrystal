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

# Import local modules
# include("<Module>.jl"); using .Module

@inline function ∇²!(∇²u, u, N, h, a)

    @tturbo for j = 2+a:N+5-a
        for i = 2+a:N+5-a
            ∇²u[i,j] = (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])/h^2
        end
    end

    return nothing

end

export ∇²!

end
