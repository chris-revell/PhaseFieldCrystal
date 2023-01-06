using DrWatson
using CairoMakie
using DataFrames

r = "17tailT_4800X_HUI_0002_0"
results = collect_results(datadir("sims", "timeResolution", r))

runID = 6

c(ϕ) = 0.5*(1-ϕ)

maskIn = load(datadir("exp_pro", "masksCompressed", r, "$r.png"))

# fibrilThreshold = -0.5# ϕ0 + 0.1
# emptySpaceThreshold = 0.5# ϕ0 - 0.1

monomerAvailability = Float64[]

for i = 1:length(results[runID, :t])
    monomerSum = 0
    for j in eachindex(results[runID, :u][i])
        if maskIn[j] > 0.5 && results[runID,:u][i][j] > 0
            monomerSum += c(results[runID,:u][i][j])
        end
    end
    push!(monomerAvailability, monomerSum)
end

fig = CairoMakie.Figure(resolution=(500,500),fontsize=32)
ax = CairoMakie.Axis(fig[1, 1], xscale=Makie.log10)
lines!(results[runID, :t][2:end], (monomerAvailability[2:end])./maximum(monomerAvailability[2:end]))
ax.xlabel = "log₁₀(Time)"
ax.ylabel = "Monomer availability"
save(datadir("sims", "timeResolution", r, "monomerAvailability$runID.png"), fig)
display(fig)