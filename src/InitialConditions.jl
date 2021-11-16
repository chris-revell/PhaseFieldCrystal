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

@views function initialConditions(imageMask,lSpace,nX,nY,ϕ0,λ,randomOrNot)

    h = lSpace/nX    # Spatial separation of grid points

    if randomOrNot == 1
        # Gaussian random field for initial u0 field
        # Lengthscale of gaussian noise (λ) set to equal lengthscale of PFC (q:=1.0)
        mean = fill(ϕ0, (nY,nX))
        cov = CovarianceFunction(2,Gaussian(λ,σ=0.1))
        ptsX = range(0, stop=lSpace, length=nX)
        ptsY = range(0, stop=h*nY, length=nY)
        grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsY, ptsX)
        # Set initial order parameter field from sample of Gaussian random field
        u0Tmp = sample(grf)
        u0Tmp .*= imageMask

    else
        u0Tmp = ones(nY,nX).*ϕ0
        u0Tmp[nX÷2-3:nX÷2+3,nY÷2-3:nY÷2+3] .+= 0.01
        u0Tmp[nX÷2-2:nX÷2+2,nY÷2-2:nY÷2+2] .+= 0.02
        u0Tmp[nX÷2-1:nX÷2+1,nY÷2-1:nY÷2+1] .+= 0.03
    end
    u0 = reshape(u0Tmp,nX*nY)

    # Allocate additional arrays for later calculations
    mat1  = zeros(nX*nY)
    mat2  = zeros(nX*nY)

return u0,mat1,mat2,h

end

export initialConditions

end
