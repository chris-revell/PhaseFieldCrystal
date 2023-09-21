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

@views function initialConditions(imageMask,lX,nX,nY,ϕ0,λ,m)

    h = lX/nX    # Spatial separation of grid points in units of q

    # Gaussian random field for initial condition 
    mean = fill(0.0, (2*nY,nX))
    cov = CovarianceFunction(2,Gaussian(λ,σ=1.0)) 
    ptsX = range(0, stop=lX, length=nX)
    ptsY = range(0, stop=2*h*nY, length=2*nY)
    grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsY, ptsX) # Gaussian random field with mean 0.0, unit variance, correlation length λ
   
    # Set initial order parameter field from sample of Gaussian random field
    u0Tmp = sample(grf)         # Form initial field by sampling Gaussian random field 
    u0Tmp1 = u0Tmp[1:nY,:] .* imageMask         # Set u0 to 0.0 for inter-cellular region
    u0Tmp2 = u0Tmp[1+nY:end,:] .* imageMask         # Set u0 to 0.0 for inter-cellular region
    u0Tmp1 .*= m                 # Multiply field by amplitude m
    u0Tmp2 .*= m                 # Multiply field by amplitude m
    u0Tmp1Mean = sum(u0Tmp1)/length(u0Tmp1[imageMask.!=0]) # Calculate actual mean of this field
    u0Tmp2Mean = sum(u0Tmp2)/length(u0Tmp2[imageMask.!=0]) # Calculate actual mean of this field
    u0Tmp1 .+= (ϕ0-u0Tmp1Mean).*imageMask # Offset field to ensure mean of ϕ0
    u0Tmp2 .+= (ϕ0-u0Tmp2Mean).*imageMask # Offset field to ensure mean of ϕ0
    u0 = zeros(Float64,2*nX*nY)
    u0[1:nX*nY] = reshape(u0Tmp1,nX*nY)   # Flatten matrix to vector 
    u0[1+nX*nY:end] = reshape(u0Tmp2,nX*nY)   # Flatten matrix to vector 

    # Pre-allocate additional arrays for use in later calculations
    mat1  = zeros(2*nX*nY)
    mat2  = zeros(2*nX*nY)

return u0,mat1,mat2,h

end

export initialConditions

end
