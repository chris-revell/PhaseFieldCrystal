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

# Import local modules
include("BoundaryConditions.jl"); using .BoundaryConditions

# Defining model as a split ode problem as per the following two links
# https://diffeq.sciml.ai/stable/solvers/split_ode_solve/
# https://diffeq.sciml.ai/stable/types/split_ode_types/#Constructors

# From Glasner, Orizaga 2016 Equation 23
# linearOperator = ((1-r+a)∇² + ∇⁶)
# f2 = ∇²(u³ - au + 2∇²u)

# f2 = ∇²(u³ - au + 2∇²u)
function f2!(du, u, p, t)

    # Unpack parameter list
    ∇², mat1, mat2, mat3, N, h, r, a = p

    # Find 2nd derivative of u
    #mat1 .= ∇²*u

    # Calculate inner component (u³ - au + 2∇²u)
    mat2 .= u.^3 .- a.*u .+ 2.0.*∇²*u

    # Find 2nd derivative of (u³ - au + 2∇²u)
    mul!(du,∇²,mat2)

    # Set values of ghost points to ensure zero flux at boundary
    #boundaryConditions!(u,N)

    return du

end

export f2!

end
