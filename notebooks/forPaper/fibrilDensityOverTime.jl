using DrWatson
using CairoMakie
using DataFrames

r = "17tailT_4800X_HUI_0002_0"
results = collect_results(datadir("sims", "timeResolution", r))

# h = 1.0
# nY = results[1, :nY]
# nX = results[1, :nX]
# ptsX = range(0, stop=h * nX, length=nX)
# ptsY = range(0, stop=h * nY, length=nY)
# for (i, t) in enumerate(results[1, :t])
#     uMat = reshape(subsetResults[1, :u][i], (nY, nX))
#     ϕSum = integrate((ptsY, ptsX), uMat)
#     aSum = integrate((ptsY, ptsX), Float64.(maskIn))
#     # push!(points,Point2(t,ϕSum/aSum))
#     display(ϕSum / aSum)
# end

maskIn = load(datadir("exp_pro", "masksCompressed", r, "$r.png"))

fibrilThreshold = -0.5# ϕ0 + 0.1
emptySpaceThreshold = 0.5# ϕ0 - 0.1

fibrilDensity = Float64[]

for i = 1:length(results[1, :t])
    monomerCount = 0
    for j in eachindex(results[1, :u][i])
        if maskIn[j] > 0.5
            if fibrilThreshold < results[1,:u][i][j] < emptySpaceThreshold
                monomerCount += 1
            end
        end
    end
    push!(fibrilDensity, monomerCount)
end

fig = CairoMakie.Figure(resolution=(500,500),fontsize=32)
ax = CairoMakie.Axis(fig[1, 1], xscale=Makie.log10)
lines!(results[1, :t][2:end], (fibrilDensity[2:end])./maximum(fibrilDensity[2:end]))
ax.xlabel = "log₁₀(Time)"
ax.ylabel = "Monomer availability"
save(datadir("sims", "timeResolution", r, "monomerAvailability.png"), fig)
display(fig)