using XLSX
using CairoMakie
using Random 
using Colors

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

xf = XLSX.readxlsx("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/Values for heatmap (collagen biosynthesis).xlsx")

t = [12.5,13.0, 13.5, 14.0, 14.5]
proteins = reverse(xf[1]["A2:A12"][:,1])
heatmapData = xf[1]["B2:F12"]

fig = CairoMakie.Figure(resolution=(600,750), fontsize=24)
ax = Axis(fig[1,1], xticks=(1:5, string.(t)), yticks=(1:11, proteins), aspect=DataAspect())
ax.xlabel = "Day"
ax.ylabel = "Protein"
heatmap!(ax,rotr90(heatmapData), colormap=:viridis)

Colorbar(fig[1,2], limits=(minimum(heatmapData), maximum(heatmapData)), colormap=:viridis)

save("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/heatmap.png",fig)

