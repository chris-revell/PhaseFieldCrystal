using DataFrames
using CairoMakie

res = collect_results(datadir("fromCSF","qPhaseSpace"))

fig = Figure(resolution=(6000,1000))

axesDict = Dict()

for (i,r) in enumerate(eachrow(subset(res, :r => r -> r .== 0.55)))
    ax = Axis(fig[1,i],aspect=DataAspect())
    uMat = reshape(r[:u][end],(r[:nY],r[:nX]))
    heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
    hidedecorations!(ax)
    hidespines!(ax)
    # Map parameters to axis in axes dictionary 
    axesDict[(r[:r],r[:ϕ0])] = ax
    Label(fig[1,i,Bottom()],"λ=$(r[:λ])",textsize=64)
end

resize_to_layout!(fig)

display(fig)