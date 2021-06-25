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

@inline @views function grad!(graduᵢ, graduⱼ, u, N, h)

    for i=1:N+5
        for j=1:N+6
            graduᵢ[i,j] = (u[i+1,j]-u[i,j])/h
            graduⱼ[j,i] = (u[j,i+1]-u[j,i])/h
        end
    end

    return nothing

end

export grad!

end
