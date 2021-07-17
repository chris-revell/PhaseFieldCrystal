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
using Laplacian

@inline function freeEnergy(sol, N, L, q, r, h)

    freeEnergies = zeros(size(sol.u))
    part1 = zeros(N+6,N+6)
    part2 = zeros(N+6,N+6)

    for (i,u) in enumerate(sol.u)

        # Find del squared u
        ∇²!(part1, u, N, h, 0)

        # Find del ^4 of u
        ∇²!(part2, part1, N, h, 1)

        # Operate in place on part2 to include additional free energy integrand terms
        part2 .+= r .+ q^4.0 .+ 2.0*q^2.0.*part1
        part2 .*= 0.5.*u
        part2 .+= 0.25.*u.^4

        # Integrate free energy matrix
        ptsX = range(0, stop=L, length=N)
        ptsY = range(0, stop=L, length=N)
        freeEnergyVal = integrate((ptsX,ptsY),part2[4:N+3,4:N+3])

        # Store free energ value from integration
        freeEnergies[i] = freeEnergyVal

    end

    display(plot(sol.t,freeEnergies))

    return nothing

end

export freeEnergy

end
