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

@inline @views function div!(divu, graduᵢ, graduⱼ, N, h)

    for i=2:N+5
        for j=2:N+5
            divu[i,j] = (graduᵢ[i,j] - graduᵢ[i-1,j] + graduⱼ[i,j] - graduⱼ[i,j-1])/h
        end
    end

    return nothing

end

export div!

end
