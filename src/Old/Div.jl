#
#  Div.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 11/06/2021.
#
#
# In place Div function

module Div

# Import Julia packages
using LinearAlgebra
using LoopVectorization

# Import local modules
# include("<Module>.jl"); using .Module

@inline function div!(divu, graduᵢ, graduⱼ, N, h)

    @tturbo for j=2:N+5
        for i=2:N+5
            divu[i,j] = (graduᵢ[i,j] - graduᵢ[i-1,j] + graduⱼ[i,j] - graduⱼ[i,j-1])/h
        end
    end

    return nothing

end

export div!

end