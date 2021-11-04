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

@views function initialConditions(lSpace,nGrid,ϕ0,λ,randomOrNot)

    if randomOrNot == 1
        # Gaussian random field for initial u0 field
        # Lengthscale of gaussian noise (λ) set to equal lengthscale of PFC (q:=1.0)
        mean = fill(ϕ0, (nGrid,nGrid))
        cov = CovarianceFunction(2,Gaussian(λ,σ=0.1))
        ptsX = range(0, stop=lSpace, length=nGrid)
        ptsY = range(0, stop=lSpace, length=nGrid)
        grf = GaussianRandomField(mean, cov, CirculantEmbedding(), ptsX, ptsY)
        # Set initial order parameter field from sample of Gaussian random field
        u0 = reshape(sample(grf),nGrid^2)
    else
        u0Tmp = ones(nGrid,nGrid).*ϕ0
        u0Tmp[nGrid÷2-3:nGrid÷2+3,nGrid÷2-3:nGrid÷2+3] .+= 0.01
        u0Tmp[nGrid÷2-2:nGrid÷2+2,nGrid÷2-2:nGrid÷2+2] .+= 0.02
        u0Tmp[nGrid÷2-1:nGrid÷2+1,nGrid÷2-1:nGrid÷2+1] .+= 0.03
        #u0Tmp[nGrid÷2,nGrid÷2] += 0.03
        u0 = reshape(u0Tmp,nGrid^2)
    end

    αVec = zeros(2*nGrid^2)
    αᵢTmp = ones(nGrid,nGrid)
    αᵢTmp[20:50,20:50] .= 0.0
    αVec[1:nGrid^2] .= reshape(αᵢTmp,nGrid^2)
    αⱼTmp = ones(nGrid,nGrid)
    αⱼTmp[20:50,20:50] .= 0.0
    αVec[1+nGrid^2:end] .= reshape(αⱼTmp,nGrid^2)
    α = spdiagm(αVec)

    # Allocate additional arrays for later calculations
    mat1  = zeros(nGrid^2)
    mat2  = zeros(nGrid^2)

    h = lSpace/nGrid    # Spatial separation of grid points

return u0,mat1,mat2,h,α

end

export initialConditions

end
