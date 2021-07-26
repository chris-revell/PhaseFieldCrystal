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
using Images

@inline @views function initialConditions(imageMask,L,N,α₀,ϕ₀,λ)

    # Gaussian random field for initial u0 field
    mean = fill(ϕ₀, (N,N))
    cov = CovarianceFunction(2,Gaussian(λ,σ=0.1))
    ptsX = range(0, stop=L, length=N)
    ptsY = range(0, stop=L, length=N)
    grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsX, ptsY)

    # Set initial order parameter field from sample of Gaussian random field
    u0 = zeros(N+6,N+6)
    u0[4:N+3,4:N+3] .= sample(grf)

    # Set spatially varying diffusivity from imported image
    αᵢ = α₀.*ones(N+5,N+6)
    αⱼ = α₀.*ones(N+6,N+5)
    for j=1:N
        for i=1:N-1
            if imageMask[i,j] == 0.0 || imageMask[i+1,j] == 0.0
                αᵢ[i+3,j+3] = 0.0
            end
        end
    end
    for j=1:N-1
        for i=1:N
            if imageMask[i,j] == 0.0 || imageMask[i,j+1] == 0.0
                αⱼ[i+3,j+3] = 0.0
            end
        end
    end

    # Allocate additional arrays for later calculations
    deriv  = zeros(N+6,N+6)
    part1  = zeros(N+6,N+6)
    part2  = zeros(N+6,N+6)
    graduᵢ = zeros(N+5,N+6)
    graduⱼ = zeros(N+6,N+5)

return u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ

end

export initialConditions

end
