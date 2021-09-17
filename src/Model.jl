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
using DifferentialEquations

# Import local modules
include("BoundaryConditions.jl"); using .BoundaryConditions
include("Laplacian.jl"); using .Laplacian
include("Grad.jl"); using .Grad
include("Div.jl"); using .Div

# Defining model as a split ode problem as per the following two links
# https://diffeq.sciml.ai/stable/solvers/split_ode_solve/
# https://diffeq.sciml.ai/stable/types/split_ode_types/#Constructors

# f1 = ((1-r+a)∇² + ∇⁶)ϕ
# f2 = ∇²(u³ - au + 2∇²u)

# From Glasner, Orizaga 2016 Equation 23

function f1(N,r,a,h)

    adj = spzeros((N+6)^2,(N+6)^2)
    dx = [- 1, 0, 1, - 1, 1, - 1, 0, 1]
    dy = [- 1, - 1, - 1, 0, 0, 1, 1, 1]
    for x=1:N+6
        for y=1:N+6
            index = (x-1)*(N+6)+y
            for ne=1:length(dx)
                newx = x + dx[ne]
                newy = y + dy[ne]
                if newx>0 && newx <=(N+6) && newy>0 && newy<=(N+6)
                    index2 = (newx-1)*(N+6) + newy
                    adj[index,index2] = 1
                end
            end
        end
    end

    degree = spdiagm(0=>sum(adj, dims=2)[:,1])

    L = (degree-adj)./h^2


    # Find (1-r+a)∇²
    linearOperator = DiffEqArrayOperator((1.0-r+a).*L .* L*L*L)

    return linearOperator,L

end

# f2 = ∇²(u³ - au + 2∇²u)
function f2!(du, u, p, t)

    # Unpack parameter list
    laplacianMatrix, mat1, mat2, mat3, N, h, r, a = p

    # Find 2nd derivative of u
    @tturbo mat1 .= laplacianMatrix*u

    # Calculate inner component (u³ - au + 2∇²u)
    @tturbo mat2 .= u.^3 .- a.*u .+ 2.0.*mat1

    # Find 2nd derivative of (u³ - au + 2∇²u)
    @tturbo mul!(du,laplacianMatrix,mat2)

    # Set values of ghost points to ensure zero flux at boundary
    #boundaryConditions!(u,N)

    return du

end

export f1, f2!

end
