using XLSX
using CairoMakie
using Random 
using Colors

xf = XLSX.readxlsx("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/Values for heatmap (collagen biosynthesis).xlsx")

t = [12.5,13.0, 13.5, 14.0, 14.5]
proteins = xf[1]["A2:A12"][:,1]
heatmapData = xf[1]["B2:F12"]

fig = CairoMakie.Figure(resolution=(600,750), fontsize=24)
ax = Axis(fig[1,1], xticks=(1:5, string.(t)), yticks=(1:11, reverse(proteins)), aspect=DataAspect())
ax.xlabel = "Day"
ax.ylabel = "Protein"
heatmap!(ax,rotr90(heatmapData), colormap=:viridis)

Colorbar(fig[1,2], limits=(minimum(heatmapData), maximum(heatmapData)), colormap=:viridis)

save("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/heatmap.png",fig)


fig2 = CairoMakie.Figure(resolution=(600,500), fontsize=24)
ax2 = Axis(fig2[1,1])#, xticks=(1:5, string.(t)))
ax2.xlabel = "Day"
ax2.ylabel = "Units?"
for i=1:11
    lines!(ax2,Point2.(t,heatmapData[i,:]),label=proteins[i],linewidth=2)
end
legend = Legend(fig2[1,2],ax2)#,tellheight = false, tellwidth = false, margin = (10, 10, 10, 10), halign = :right, valign = :top)#, orientation = :horizontal

save("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/heatmapAsLines.png",fig2)

