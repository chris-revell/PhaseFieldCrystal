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
using SparseArrays
using LinearAlgebra

@views function initialConditions(imageMask,lSpace,nGrid,ϕ0,λ,randomOrNot)

    if randomOrNot == 1
        # Gaussian random field for initial u0 field
        # Lengthscale of gaussian noise (λ) set to equal lengthscale of PFC (q:=1.0)
        mean = fill(ϕ0, (nGrid,nGrid))
        cov = CovarianceFunction(2,Gaussian(λ,σ=0.1))
        ptsX = range(0, stop=lSpace, length=nGrid)
        ptsY = range(0, stop=lSpace, length=nGrid)
        grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsX, ptsY)
        # Set initial order parameter field from sample of Gaussian random field
        u0Tmp = sample(grf)
        u0Tmp .*= imageMask

    else
        u0Tmp = ones(nGrid,nGrid).*ϕ0
        u0Tmp[nGrid÷2-3:nGrid÷2+3,nGrid÷2-3:nGrid÷2+3] .+= 0.01
        u0Tmp[nGrid÷2-2:nGrid÷2+2,nGrid÷2-2:nGrid÷2+2] .+= 0.02
        u0Tmp[nGrid÷2-1:nGrid÷2+1,nGrid÷2-1:nGrid÷2+1] .+= 0.03
    end
    u0 = reshape(u0Tmp,nGrid^2)

    # Allocate additional arrays for later calculations
    mat1  = zeros(nGrid^2)
    mat2  = zeros(nGrid^2)

    h = lSpace/nGrid    # Spatial separation of grid points

return u0,mat1,mat2,h

end

export initialConditions

end
