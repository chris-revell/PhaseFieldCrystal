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

@inline function freeEnergy(sol, laplacianMatrix, mat1, mat2, N, L, r)

    freeEnergies = zeros(size(sol.u))
    mat1 = zeros((N+6)^2)
    mat2 = zeros((N+6)^2)
    mat3 = zeros(N+6,N+6)

    for (i,u) in enumerate(sol.u)

        # Find del ^4 of u
        mat1 .= laplacianMatrix*u
        mat2 .= laplacianMatrix*mat1

        # Operate in place on part2 to include additional free energy integrand terms
        mat2 .+= r .+ 2.0.*mat1
        mat2 .*= 0.5.*u
        mat2 .+= 0.25.*u.^4

        # Integrate free energy matrix
        ptsX = range(0, stop=L, length=N)
        ptsY = range(0, stop=L, length=N)
        mat3 .= reshape(mat2,(N+6,N+6))
        freeEnergyVal = integrate((ptsX,ptsY),mat3[4:N+3,4:N+3])

        # Store free energ value from integration
        freeEnergies[i] = freeEnergyVal

    end

    return freeEnergies

end

export freeEnergy

end
