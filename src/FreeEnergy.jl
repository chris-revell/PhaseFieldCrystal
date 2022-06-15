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

@inline function freeEnergy(sol, ∇², mat1, mat2, nX, nY, h, r)

    freeEnergies = zeros(size(sol.u))

    for (i,u) in enumerate(sol.u)

        # Find del ^4 of u
        mat1 .= ∇²*u
        mat2 .= ∇²*mat1

        # Operate in place on part2 to include additional free energy integrand terms
        mat2 .+= r .+ 2.0.*mat1
        mat2 .*= 0.5.*u
        mat2 .+= 0.25.*u.^4

        # Integrate free energy matrix
        ptsX = range(0, stop=h*nX, length=nX)
        ptsY = range(0, stop=h*nY, length=nY)
        mat3 = reshape(mat2,(nY,nX))
        freeEnergyVal = integrate((ptsY,ptsX),mat3)

        # Store free energy value from integration
        freeEnergies[i] = freeEnergyVal

    end

    return freeEnergies

end

export freeEnergy

end
