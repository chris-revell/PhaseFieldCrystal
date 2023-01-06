using DrWatson;
@quickactivate;
using DataFrames
using CairoMakie
using Colors
using FromFile
using GeometryBasics
using NumericalIntegration
using Printf

r = "17tailT_4800X_HUI_0002_0"
runID = 6

results = collect_results(datadir("sims", "timeResolution", r))

c(ϕ) = 0.5 * (1 - ϕ)

fig = CairoMakie.Figure(figure_padding=30,resolution=(500, 500), fontsize=32)
ax2 = CairoMakie.Axis(fig[1, 1], xscale=Makie.log10, yticklabelcolor=:blue, ylabelcolor=:blue)
lines!(ax2, results[runID, :t][2:end], results[runID, :freeEnergies][2:end] ./ 1E7, color=(:blue, 0.75))
ax2.xlabel = "Time"
ax2.ylabel = "Free Energy /10⁷"


maskIn = load(datadir("exp_pro", "masksCompressed", r, "$r.png"))
# fibrilThreshold = -0.5# ϕ0 + 0.1
# emptySpaceThreshold = 0.5# ϕ0 - 0.1
monomerAvailability = Float64[]
for i = 1:length(results[runID, :t])
    monomerSum = 0
    for j in eachindex(results[runID, :u][i])
        if maskIn[j] > 0.5 && results[runID, :u][i][j] > 0
            monomerSum += c(results[runID, :u][i][j])
        end
    end
    push!(monomerAvailability, monomerSum)
end
ax = CairoMakie.Axis(fig[1, 1], xscale=Makie.log10, yticklabelcolor=:green, ylabelcolor=:green)
lines!(results[runID, :t][2:end], (monomerAvailability[2:end]) ./ maximum(monomerAvailability[2:end]), color=(:green, 0.75))
ax.xlabel = "log₁₀(Time)"
ax.ylabel = "Monomer availability"

ax.yaxisposition = :right
ax.yticklabelalign = (:left, :center)
ax.xticklabelsvisible = false
ax.xticklabelsvisible = false
ax.xlabelvisible = false

linkxaxes!(ax, ax2)



colsize!(fig.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig)
display(fig)
save(datadir("sims", "timeResolution", r, "freeEnergyAndMonomerAvailabilityLog$runID.png"), fig)