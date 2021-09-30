#
#  Model.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#

module Model

# Import Julia packages
using LinearAlgebra
using LoopVectorization
using SparseArrays
using Octavian

# Import local modules
include("BoundaryConditions.jl"); using .BoundaryConditions

# Defining model as a split ode problem as per the following two links
# https://diffeq.sciml.ai/stable/solvers/split_ode_solve/
# https://diffeq.sciml.ai/stable/types/split_ode_types/#Constructors

# From Glasner, Orizaga 2016 Equation 23
# linearOperator = ((1-r+a)∇² + ∇⁶)
# f2 = ∇²(u³ - au + 2∇²u)

function f2!(du, u, p, t)

    # Unpack parameter list
    ∇², linearOperator, mat1, mat2, nGrid, h, r, a = p

    # Find 2nd derivative of u
    matmul!(mat1,∇²,u)

    # Calculate inner component (u³ - au + 2∇²u)
    @tturbo mat2 .= u.^3 .- a.*u .+ 2.0.*mat1

    # Find 2nd derivative of (u³ - au + 2∇²u)
    matmul!(du,∇²,mat2)

    matmul!(mat2,linearOperator,u)

    du .+= mat2

    # Set values of ghost points to ensure zero flux at boundary
    boundaryConditions!(u, nGrid, 3)

    duTmp = reshape(du,(nGrid+6,nGrid+6))
    duTmp[:,1] .= 0.0
    duTmp[:,2] .= 0.0
    duTmp[:,3] .= 0.0
    duTmp[1,:] .= 0.0
    duTmp[2,:] .= 0.0
    duTmp[3,:] .= 0.0
    duTmp[:,nGrid+3+1] .= 0.0
    duTmp[:,nGrid+3+2] .= 0.0
    duTmp[:,nGrid+3+3] .= 0.0
    duTmp[nGrid+3+1,:] .= 0.0
    duTmp[nGrid+3+2,:] .= 0.0
    duTmp[nGrid+3+3,:] .= 0.0

    return du

end

function PFC!(du, u, p, t)

    # Unpack parameter list
    ∇², mat1, mat2, r = p
    
    # Find Laplacian of u
    mat1 = ∇²*u

    # Calculate inner component (∇²ϕ + q²ϕ)
    mat1 .+= u

    # Find Laplacian of (∇²ϕ + q²ϕ)
    mat2 .= ∇²*mat1

    # Calculate full term within outermost Laplacian (rϕ + ∇²(∇²ϕ + q²ϕ) + q⁴ + ϕ³) = rϕ + ∇²(mat1) + q⁴ + ϕ³
    mat2 .+= r.*u .+ u.^3

    du .= ∇²*mat2

    return du

end

function diffusion!(du, u, p, t)

    # Unpack parameter list
    ∇², mat1, mat2, r = p

    du .= -1.0.*∇²*u

    return du

end

function cahnHilliard!(du, u, p, t)

    # Unpack parameter list
    ∇², mat1, mat2, r = p

    mat1 .= -0.00005.*∇²*u 
    mat1.+= u.^3 .- u

    du .= ∇²*mat1

    return du

end

export f2!, PFC!, diffusion!, cahnHilliard!

end
