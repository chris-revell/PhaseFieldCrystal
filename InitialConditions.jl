#
#  InitialConditions.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 13/05/2021.
#
#
# Function to define necessary matrices and set initial conditions for order parameter and mobility fields

module InitialConditions

using GaussianRandomFields
using NumericalIntegration

@inline @views function initialConditions(L,N,α₀,ϕ₀,λ)

    # Gaussian random field for initial u0 field
    mean = fill(ϕ₀, (N,N))
    cov = CovarianceFunction(2,Gaussian(λ,σ=0.1))
    ptsX = range(0, stop=L, length=N)
    ptsY = range(0, stop=L, length=N)
    grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsX, ptsY)

    # Set initial order parameter field from sample of Gaussian random field
    u0 = zeros(N+6,N+6)
    u0[4:N+3,4:N+3] .= sample(grf)

    # Integrate over domain to find actual average order parameter (ϕ₀ mean in random distribution, but randomness means exact mean may differ)
    ϕ₀Real = integrate((ptsX,ptsY),u0[4:N+3,4:N+3]/(L*L))

    # Set spatially varying diffusivity
    αᵢ = α₀.*ones(N+5,N+6)
    αⱼ = α₀.*ones(N+6,N+5)

    # Want to block off i = 11:20, j=101:150
    ilow = 35
    ihigh = 70
    jlow = 30
    jhigh = 80
    αᵢ[ilow-1:ihigh,jlow:jhigh] .= 0.0
    αⱼ[ilow:ihigh,jlow-1:jhigh] .= 0.0

    # # Set smooth edges
    # αᵢ[ilow-2,jlow-1:jhigh+1] .= 0.1
    # αᵢ[ihigh+1,jlow-1:jhigh+1] .= 0.1
    # αᵢ[ilow-1:ihigh+1,jlow-2] .= 0.1
    # αᵢ[ilow-1:ihigh+1,jhigh+1] .= 0.1
    # αᵢ[ilow-3,jlow-2:jhigh+2] .= 0.2
    # αᵢ[ihigh+2,jlow-2:jhigh+2] .= 0.2
    # αᵢ[ilow-2:ihigh+2,jlow-3] .= 0.2
    # αᵢ[ilow-2:ihigh+2,jhigh+2] .= 0.2
    # αᵢ[ilow-4,jlow-3:jhigh+3] .= 0.3
    # αᵢ[ihigh+3,jlow-3:jhigh+3] .= 0.3
    # αᵢ[ilow-3:ihigh+3,jlow-4] .= 0.3
    # αᵢ[ilow-3:ihigh+3,jhigh+3] .= 0.3
    # αᵢ[ilow-5,jlow-4:jhigh+4] .= 0.4
    # αᵢ[ihigh+4,jlow-4:jhigh+4] .= 0.4
    # αᵢ[ilow-4:ihigh+4,jlow-5] .= 0.4
    # αᵢ[ilow-4:ihigh+4,jhigh+4] .= 0.4

    # Allocate additional arrays for later calculations
    deriv  = zeros(N+6,N+6)
    part1  = zeros(N+6,N+6)
    part2  = zeros(N+6,N+6)
    graduᵢ = zeros(N+5,N+6)
    graduⱼ = zeros(N+6,N+5)

return u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ,ϕ₀Real

end

export initialConditions

end
