using XLSX
using CairoMakie
using Random 
using Colors

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

xf = XLSX.readxlsx("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/extractedMSvalues.xlsx")

t = [12.5,13.0, 13.5, 14.0, 14.5]

colors = [:black, :red, :green, :blue, :yellow]

fig = CairoMakie.Figure(resolution=(750,750), fontsize=32)
ax = Axis(fig[1,1])
ax.xlabel = "Day"
ax.ylabel = "log₂ of fold change in expression vs E12.5"
for i=1:4
    means = [0.0]
    sems = [0.0]
    datalabel = xf[1][1,(i-1)*3+2]
    for j=3:6
        push!(means,xf[1][j,(i-1)*3+2])
        push!(sems,xf[1][j,(i-1)*3+3])
    end    
    lines!(ax,t,means,color=colors[i])
    errorbars!(ax,t,means,sems,label=datalabel,color=colors[i])
end
legend = Legend(fig[1,1],ax,tellheight = false, tellwidth = false, margin = (10, 10, 10, 10), halign = :right, valign = :top)#, orientation = :horizontal

save("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/extractedMSvalues1.png",fig)

fig = CairoMakie.Figure(resolution=(750,750), fontsize=32)
ax = Axis(fig[1,1])
ax.xlabel = "Day"
ax.ylabel = "log₂ of fold change in expression vs E12.5"
for i=5:7
    means = [0.0]
    sems = [0.0]
    datalabel = xf[1][1,(i-1)*3+2]
    for j=3:6
        push!(means,xf[1][j,(i-1)*3+2])
        push!(sems,xf[1][j,(i-1)*3+3])
    end    
    lines!(ax,t,means,color=colors[i-4])
    errorbars!(ax,t,means,sems,label=datalabel,color=colors[i-4])
end
legend = Legend(fig[1,1],ax,tellheight = false, tellwidth = false, margin = (10, 10, 10, 10), halign = :right, valign = :top)#, orientation = :horizontal
# display(fig)
save("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/extractedMSvalues2.png",fig)