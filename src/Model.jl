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
using SparseArrays

function splitNonlinearPart!(du, u, p, t)
    # Defining model as a split ode problem as per the following two links
    # https://diffeq.sciml.ai/stable/solvers/split_ode_solve/
    # https://diffeq.sciml.ai/stable/types/split_ode_types/#Constructors

    # From Glasner, Orizaga 2016 Equation 23
    # linearOperator = ((1-r+a)∇² + ∇⁶)
    # f2 = ∇²(u³ - au + 2∇²u)

    # Unpack parameter list
    ∇², linearOperator, mat1, mat2, ∇x, ∇y, divX, divY, r, a, αᵢ, αⱼ = p

    # Find 2nd derivative of u
    mul!(mat1,∇²,u)
    # Calculate inner component (u³ - au + 2∇²u)
    mat1 .*= 2.0
    mat1 .+= u.^3 .- a.*u

    # Find 2nd derivative of (u³ - au + 2∇²u)
    #du .= divX*(αⱼ*∇x)*mat1 .+ divY*(αᵢ*∇y)*mat1
    mul!(du,∇²,mat1)

    return du
end

function PFC!(du, u, p, t)
    # Unpack parameter list
    ∇², linearOperator, mat1, mat2, r, a = p
    # Find Laplacian of u
    mul!(mat1,∇²,u)
    # Calculate inner component (∇²ϕ + q²ϕ)
    mat1 .+= u
    # Find Laplacian of (∇²ϕ + q²ϕ)
    mul!(mat2,∇²,mat1)
    # Calculate full term within outermost Laplacian (rϕ + ∇²(∇²ϕ + q²ϕ) + q⁴ + ϕ³) = rϕ + ∇²(mat1) + q⁴ + ϕ³
    mat2 .+= u.^3 .- r.*u
    mul!(du,∇²,mat2)
    return du
end

export splitNonlinearPart!, PFC!

end
