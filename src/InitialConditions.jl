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

@views function initialConditions(imageMask,m,lX,nX,nY,ϕ0,λ)

    h = lX/nX    # Spatial separation of grid points

    # Gaussian random field for initial u0 field
    # Lengthscale of gaussian noise (λ) set to equal lengthscale of PFC (q:=1.0)
    mean = fill(ϕ0, (nY,nX))
    cov = CovarianceFunction(2,Gaussian(λ,σ=m))
    ptsX = range(0, stop=lX, length=nX)
    ptsY = range(0, stop=h*nY, length=nY)
    grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsY, ptsX)
    # Set initial order parameter field from sample of Gaussian random field
    u0Tmp = sample(grf)
    u0Tmp .*= imageMask

    u0 = reshape(u0Tmp,nX*nY)

    # Allocate additional arrays for later calculations
    mat1  = zeros(nX*nY)
    mat2  = zeros(nX*nY)

return u0,mat1,mat2,h

end

export initialConditions

end
