#
#  FreeEnergy.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 23/07/2021.
#
#
# Integrate free energy over domain
# F[ϕ] ≡ ∫dx[(ϕ/2)(r + (q² + ∇²)²)ϕ + ϕ⁴/4]
# F[ϕ] ≡ ∫dx[(ϕ/2)(r + q⁴ + q²∇²ϕ + ∇⁴ϕ) + ϕ⁴/4]

module FreeEnergy

# Import Julia packages
using LinearAlgebra
using NumericalIntegration

# Import local Julia modules
using Laplacian

@inline function freeEnergy(u, N, L, q, r, h)

    part1 = zeros(N+6,N+6)
    part2 = zeros(N+6,N+6)

    # Find del squared u
    ∇²!(part1, u, N, h, 0)

    # Find del ^4 of u
    ∇²!(part2, part1, N, h, 1)

    part2 .+= r .+ q^4.0 .+ q^2.0.*part1

    part2 .*= 0.5.*u

    part2 .+= 0.25.*u.^4

    ptsX = range(0, stop=L, length=N)
    ptsY = range(0, stop=L, length=N)
    freeEnergyVal = integrate((ptsX,ptsY),part2[4:N+3,4:N+3])

    return freeEnergyVal

end

export freeEnergy

end
