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

@views function initialConditions(imageMask,L,N,α₀,ϕ₀,λ,nGhosts)

    # Gaussian random field for initial u0 field
    mean = fill(ϕ₀, (N,N))
    cov = CovarianceFunction(2,Gaussian(λ,σ=0.1))
    ptsX = range(0, stop=L, length=N)
    ptsY = range(0, stop=L, length=N)
    grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsX, ptsY)

    # Set initial order parameter field from sample of Gaussian random field
    u0 = zeros(N+2*nGhosts,N+2*nGhosts)
    u0[nGhosts+1:N+nGhosts,nGhosts+1:N+nGhosts] .= sample(grf)

    # Set spatially varying diffusivity from imported image
    αᵢ = α₀.*ones(N+2*nGhosts-1,N+2*nGhosts)
    αⱼ = α₀.*ones(N+2*nGhosts,N+2*nGHosts-1)
    for j=1:N
        for i=1:N-1
            if imageMask[i,j] == 0.0 || imageMask[i+1,j] == 0.0
                αᵢ[i+nGhosts,j+nGhosts] = 0.0
            end
        end
    end
    for j=1:N-1
        for i=1:N
            if imageMask[i,j] == 0.0 || imageMask[i,j+1] == 0.0
                αⱼ[i+nGhosts,j+nGhosts] = 0.0
            end
        end
    end

    # Allocate additional arrays for later calculations
    deriv  = zeros(N+2*nGhosts,N+2*nGhosts)
    part1  = zeros(N+2*nGhosts,N+2*nGhosts)
    part2  = zeros(N+2*nGhosts,N+2*nGhosts)
    graduᵢ = zeros(N+2*nGhosts-1,N+2*nGhosts)
    graduⱼ = zeros(N+2*nGhosts,N+2*nGhosts-1)

return u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ

end

export initialConditions

end
