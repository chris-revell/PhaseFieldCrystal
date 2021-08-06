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
using LoopVectorization

@inline function ∇²!(∇²u, u, N, h, a, nGhosts)

    @tturbo for j=2+a:N+(nGhosts*2-1)-a
        for i=2+a:N+(nGhosts*2-1)-a
            ∇²u[i,j] = (-u[i-2,j]/12.0 + 4.0*u[i-1,j]/3.0 - 5.0*u[i,j]/2.0 + 4.0*u[i+1,j]/3.0 - u[i+2,j]/12.0 +
                            -u[i,j-2]/12.0 + 4.0*u[i,j-1]/3.0 - 5.0*u[i,j]/2.0 + 4.0*u[i,j+1]/3.0 - u[i,j+2]/12.0)/h^2
        end
    end

    return nothing

end

export ∇²!

end
