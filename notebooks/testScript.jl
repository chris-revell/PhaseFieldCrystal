
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using DrWatson
using IncidenceMatrixTools
using GaussianRandomFields
using FastBroadcast
using SciMLOperators
using CairoMakie
using Printf
using Statistics
# using NamedTuples

dims = [250,250]
spaceLengths = [100.0, 100.0]
spacing = spaceLengths./(dims.-1) 
m = 0.1
a = 2.0
c = 0.2
λ = 0.1*spaceLengths[1]
δt = 0.1
r = 0.8
ϕ₀ = 0.4
q = 1.0
m = 0.1
σ = 0.1
𝒟 = 100.0
tMax = 100.0

cov = CovarianceFunction(length(dims), Gaussian(λ, σ=σ))
pts1 = range(0, stop=spaceLengths[1], length=dims[1])
pts2 = range(0, stop=spaceLengths[2], length=dims[2])
grf = GaussianRandomField(cov, CirculantEmbedding(), pts1, pts2, minpadding=113)
u0mat = m.*sample(grf)[1:dims[1], 1:dims[2]]
u0Mean = mean(u0mat)
u0mat .+= (ϕ₀-u0Mean)
u0 = reshape(u0mat, prod(dims))

# fig1 = Figure(figure_padding=0, size=(1000, 1000), fontsize=64)
# ax1 = CairoMakie.Axis(fig1[1, 1], aspect=DataAspect())
# uInternal1 = Observable(zeros(dims...))
# heatmap!(ax1, u0mat)
# display(fig1)

# Incidence matrices
A   = makeIncidenceMatrix3D(dims, periodicity=[1,1,0])
Aᵀ  = transpose(A)
# Weights
W   = vertexVolumeWeightsMatrix(dims, spacing)
W⁻¹ =  vertexVolumeWeightsInverseMatrix(dims, spacing)
L⁻¹ = edgeLengthInverseMatrix(dims, spacing, periodicity=[1,1])
# Gradient operators
∇ₑ = L⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 
# Diffusivity field over edges 
Aperpₑ = edgePerpendicularAreaMatrix(dims, spacing, periodicity=[1,1])
𝒟ₑ     = 𝒟.*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

# Number of edges over each dimension 
# dimEdgeCount = Int64[]
# for i=1:length(dims)
#     push!(dimEdgeCount, (dims[i]-1)*prod(dims[Not(i)]))
# end
# nVerts  = prod(dims)          # Total number of vertices 
# nEdges  = sum(dimEdgeCount)   # Total number of edges over all dimensions 
# Matrices for picking out ν and xy directions in derivatives 
# Pν  = spdiagm(vcat(ones(Int64, dimEdgeCount[1]), zeros(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all xy edges 
# Pxy  = spdiagm(vcat(zeros(Int64, dimEdgeCount[1]), ones(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all ν edges 

# Create matrix for linear component of PFC equation
constantComponent = spdiagm(ones(Float64, prod(dims)) .- r .+ a )
ℳ = ∇cdot*𝒟ₑ*∇ₑ*(constantComponent .+ ∇cdot*∇ₑ*∇cdot*∇ₑ)

function splitNonlinearPart!(du, u, p, t)
    # Defining model as a split ode problem as per the following two links
    # https://diffeq.sciml.ai/stable/solvers/split_ode_solve/
    # https://diffeq.sciml.ai/stable/types/split_ode_types/#Constructors
    # From Glasner, Orizaga 2016 Equation 23
    # ℳ = ((1-r+a)∇² + ∇⁶)
    # f2 = ∇²(u³ - au + 2∇²u)
    # Unpack parameter list
    @unpack ∇², ∇α∇, ℳ, tmpVec, r, a, 𝒟ₑ = p
    # Find 2nd derivative of u
    mul!(tmpVec,∇²,u)
    # Calculate inner component (u³ - au + 2∇²u)
    @.. thread=true tmpVec .*= 2.0
    @.. thread=false tmpVec .+= u.^3 .- a.*u
    # Find 2nd derivative of (u³ - au + 2∇²u)
    mul!(du, ∇α∇, tmpVec)
    return du
end

p = (
    ∇² = ∇cdot*∇ₑ,
    ∇α∇ = ∇cdot*𝒟ₑ*∇ₑ,
    ℳ = ℳ,
    tmpVec = zeros(prod(dims)), # Pre-allocate additional arrays for use in later calculations
    r = r,
    a = a,
    𝒟ₑ = 𝒟ₑ,
)

# Define split ODE problem
prob = SplitODEProblem(MatrixOperator(ℳ), splitNonlinearPart!, reshape(u0, prod(dims)), (0.0, tMax), p)
sol = solve(prob, ETDRK2(krylov=true, m=50), dt=δt, saveat = tMax/100, progress=true)

#%%

fig1 = Figure(figure_padding=0, size=(1000, 1000), fontsize=64)
ax1 = CairoMakie.Axis(fig1[1, 1], aspect=DataAspect())
uInternal1 = Observable(zeros(dims...))
heatmap!(ax1, uInternal1, colorrange=(-1.0, 1.0), colormap=(:bwr, 1.0))
hidedecorations!(ax1)
hidespines!(ax1)
ax1.title = "t=0.0"
ax1.yreversed = true
resize_to_layout!(fig1)
mov = VideoStream(fig1, framerate=10)
for i=1:length(sol.t)
    display(i)
    ax1.title = "t=$(@sprintf("%.2f", sol.t[i]))"
    uInternal1[] .= transpose(reshape(sol.u[i], (dims...)))
    uInternal1[] = uInternal1[]
    recordframe!(mov)
end
safesave(datadir("sims", "tstNew.mp4"), mov)


fig1 = Figure(figure_padding=0, size=(1000, 1000), fontsize=64)
ax1 = CairoMakie.Axis(fig1[1, 1], aspect=DataAspect())
uInternal1 = Observable(zeros(dims...))
heatmap!(ax1, u0mat, colorrange=(-1,1), colormap=(:bwr,1.0))

Colorbar(fig1[1,2], colorrange=(-1,1), colormap=(:bwr,1.0))
display(fig1)