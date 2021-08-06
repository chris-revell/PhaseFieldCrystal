#
#  Div.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 11/06/2021.
#
#
# In place Div function

module Div

using LinearAlgebra
using LoopVectorization

@inline function div!(divu, graduᵢ, graduⱼ, N, h, nGhosts)

    @tturbo for j=2:N+nGhosts*2-1
        for i=2:N+nGhosts*2-1
            divu[i,j] = (graduᵢ[i,j] - graduᵢ[i-1,j] + graduⱼ[i,j] - graduⱼ[i,j-1])/h
        end
    end

    return nothing

end

export div!

end
