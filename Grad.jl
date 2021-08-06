#
#  Grad.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 11/06/2021.
#
#
# In place Grad function

module Grad

using LinearAlgebra
using LoopVectorization

@inline function grad!(graduᵢ, graduⱼ, u, αᵢ, αⱼ, N, h, nGhosts)

    @tturbo for j=1:N+nGhosts*2
        for i=1:N+nGhosts*2-1
            graduᵢ[i,j] = αᵢ[i,j]*(u[i+1,j]-u[i,j])/h
        end
    end

    @tturbo for i=1:N+nGhosts*2-1
        for j=1:N+nGhosts*2
            graduⱼ[j,i] = αⱼ[j,i]*(u[j,i+1]-u[j,i])/h
        end
    end

    return nothing

end

export grad!

end
