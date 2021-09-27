#
#  FreeEnergy.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 23/07/2021.
#
#
# Integrate free energy over domain
# F[ϕ] ≡ ∫dx[(ϕ/2)(r + (q² + ∇²)²)ϕ + ϕ⁴/4]
# F[ϕ] ≡ ∫dx[(ϕ/2)(r + q⁴ + 2q²∇²ϕ + ∇⁴ϕ) + ϕ⁴/4]

module FreeEnergy

# Import Julia packages
using LinearAlgebra
using NumericalIntegration
using Plots

# Import local Julia modules
include("Laplacian.jl"); using .Laplacian

@inline function freeEnergy(sol, ∇², mat1, mat2, nGrid, nGhosts, lSpace, r)

    freeEnergies = zeros(size(sol.u))
    mat1 = zeros((nGrid+nGhosts)^2)
    mat2 = zeros((nGrid+nGhosts)^2)
    mat3 = zeros(nGrid+nGhosts,nGrid+nGhosts)

    for (i,u) in enumerate(sol.u)

        # Find del ^4 of u
        mat1 .= ∇²*u
        mat2 .= ∇²*mat1

        # Operate in place on part2 to include additional free energy integrand terms
        mat2 .+= r .+ 2.0.*mat1
        mat2 .*= 0.5.*u
        mat2 .+= 0.25.*u.^4

        # Integrate free energy matrix
        ptsX = range(0, stop=lSpace, length=nGrid)
        ptsY = range(0, stop=lSpace, length=nGrid)
        mat3 .= reshape(mat2,(nGrid+nGhosts,nGrid+nGhosts))
        freeEnergyVal = integrate((ptsX,ptsY),mat3[4:nGrid+3,4:nGrid+3])

        # Store free energ value from integration
        freeEnergies[i] = freeEnergyVal

    end

    return freeEnergies

end

export freeEnergy

end
