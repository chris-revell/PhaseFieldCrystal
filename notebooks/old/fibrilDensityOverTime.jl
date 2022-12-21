using DrWatson
using CairoMakie

fileName = "data/sims/timeResolution/testMask/lX=200.0m=0.1maskFileName=testMasknX=500nY=383r=0.7tMax=2000.0δt=0.5λ=1.0ϕ0=0.42.jld2"

dataIn = load(fileName)
# @unpack u, t, ϕ0, r, m, nX, nY, lX, a, δt, tMax = data
@unpack ϕ0, u, t = dataIn

maskFileName = datadir("exp_pro/testMask.png")
maskImage = load(maskFileName)#(datadir("exp_pro","masksCompressed",maskFileName,"$maskFileName.png"))
maskSpaceSize = length(filter(x->x>0.5,maskImage))

fibrilThreshold = ϕ0+0.1
emptySpaceThreshold = ϕ0-0.1

fibrilDensity = Float64[]

for i=1:length(t)
    fibrilCount = 0
    for j in eachindex(u[1])
        if maskImage[j]>0.5
            if emptySpaceThreshold<u[i][j]<fibrilThreshold
                fibrilCount+=1
            end
        end
    end
    push!(fibrilDensity,fibrilCount/maskSpaceSize)
end 

fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1],xscale = Makie.pseudolog10)
lines!(t,fibrilDensity)
ax.xlabel = L"log_{10}(Time)"
ax.ylabel = "Monomer availability"
# ax.ylabel = "Fibril density"
display(fig)
save("test.png",fig)