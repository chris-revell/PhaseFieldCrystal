#
#  InitialConditions.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 13/05/2021.
#
#
# Function to define necessary matrices and set initial conditions for order parameter and mobility fields

module InitialConditions

# Import Julia packages
using GaussianRandomFields
using NumericalIntegration
using Images

# Import local modules
include("BoundaryConditions.jl"); using .BoundaryConditions

@views function initialConditions(lSpace,nGrid,nGhosts,ϕ₀,λ)

    # Gaussian random field for initial u0 field
    # Lengthscale of gaussian noise (λ) set to equal lengthscale of PFC (q:=1.0)
    mean = fill(ϕ₀, (nGrid,nGrid))
    cov = CovarianceFunction(2,Gaussian(λ,σ=0.1))
    ptsX = range(0, stop=lSpace, length=nGrid)
    ptsY = range(0, stop=lSpace, length=nGrid)
    grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsX, ptsY)

    # Set initial order parameter field from sample of Gaussian random field
    u0 = zeros(nGrid+nGhosts,nGrid+nGhosts)
    inset = nGhosts÷2
    u0[inset+1:inset+nGrid,inset+1:inset+nGrid] .= sample(grf)

    # Set values of ghost points to ensure zero flux at boundary
    boundaryConditions!(u0,nGrid,nGhosts÷2)

    u0 = reshape(u0,(nGrid+nGhosts)^2)

    # Allocate additional arrays for later calculations
    mat1  = zeros((nGrid+nGhosts)^2)
    mat2  = zeros((nGrid+nGhosts)^2)
    mat3  = zeros((nGrid+nGhosts)^2)

    h = lSpace/(nGrid-1)    # Spatial separation of grid points

return u0,mat1,mat2,mat3,h

end

export initialConditions

end
