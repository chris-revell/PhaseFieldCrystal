using DrWatson
using CairoMakie

r = "17tailT_4800X_HUI_0002_0"
voronoiSizeThresh = 1.3


h = 1.0#subsetResults[1,:h]
nY = subsetResults[1, :nY]
nX = subsetResults[1, :nX]
ptsX = range(0, stop=h * nX, length=nX)
ptsY = range(0, stop=h * nY, length=nY)
for (i, t) in enumerate(subsetResults[1, :t])
    uMat = reshape(subsetResults[1, :u][i], (nY, nX))
    ϕSum = integrate((ptsY, ptsX), uMat)
    aSum = integrate((ptsY, ptsX), Float64.(maskIn))
    # push!(points,Point2(t,ϕSum/aSum))
    display(ϕSum / aSum)
end



results = collect_results(datadir("sims", "timeResolution", r))
nY = results[1, :nY]
nX = results[1, :nX]

maskIn = load(datadir("exp_pro", "masksCompressed", r, "$r.png"))
maskImage = fill(RGBA(1, 1, 1, 1), size(maskIn))
for i = 1:size(maskIn)[1]
    for j = 1:size(maskIn)[2]
        if maskIn[i, j] < 0.5
            maskImage[i, j] = RGBA(0.0, 0.0, 0.0, 1.0)
        else
            maskImage[i, j] = RGBA(0.0, 0.0, 0.0, 0.0)
        end
    end
end
maskSpaceSize = length(filter(x -> x > 0.5, maskImage))

fibrilThreshold = -0.5# ϕ0 + 0.1
emptySpaceThreshold = 0.5# ϕ0 - 0.1

fibrilDensity = Float64[]

for i = 1:length(t)
    fibrilCount = 0
    for j in eachindex(results[1, :u][i])
        if maskImage[j] > 0.5
            if emptySpaceThreshold < u[i][j] < fibrilThreshold
                fibrilCount += 1
            end
        end
    end
    push!(fibrilDensity, fibrilCount / maskSpaceSize)
end

fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1, 1], xscale=Makie.pseudolog10)
lines!(t, fibrilDensity)
ax.xlabel = L"log_{10}(Time)"
ax.ylabel = "Monomer availability"
# ax.ylabel = "Fibril density"
display(fig)
save("test.png", fig)