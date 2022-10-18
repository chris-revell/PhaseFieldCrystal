using DrWatson
using CairoMakie

fileName = "data/fromCSF/longRuns/18tailT_4800X_HUI_0003_0/lX=246.0m=0.1nX=614nY=456r=0.7tMax=10000.0δt=0.1λ=2.0ϕ0=0.41.jld2"

dataIn = load(fileName)
# @unpack u, t, ϕ0, r, m, nX, nY, lX, a, δt, tMax = data
@unpack ϕ0, u, t = dataIn

maskFileName = "18tailT_4800X_HUI_0003_0"
maskImage = load(datadir("exp_pro","masksCompressed",maskFileName,"$maskFileName.png"))
maskSpaceSize = length(filter(x->x>0.5,maskImage))

fibrilThreshold = ϕ0

fibrilDensity = Float64[]

for i=1:length(t)
    fibrilCount = 0
    for j in eachindex(u[1])
        if maskImage[j]>0.5
            if u[i][j]>fibrilThreshold
                fibrilCount+=1
            end
        end
    end
    push!(fibrilDensity,fibrilCount/maskSpaceSize)
end 

fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1])
lines!(t,fibrilDensity)
ax.xlabel = "Time"
ax.ylabel = "Fibril density"
display(fig)
