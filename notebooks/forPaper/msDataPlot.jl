using XLSX
using CairoMakie
using Random 
using Colors

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

xf = XLSX.readxlsx("/Users/christopher/Dropbox (The University of Manchester)/RevellHerrera1/MS_Data/extractedMSvalues.xlsx")

fig = CairoMakie.Figure(resolution=(800,500))
ax = Axis(fig[1,1])
ax.xlabel = "Day /E"
ax.ylabel = "Percentage change"

t = [13.0, 13.5, 14.0, 14.5]#(xf[1]["A3:A6"])[:,1]

for i=1:7
    means = Float64[]
    sems = Float64[]
    datalabel = xf[1][1,(i-1)*3+2]
    for j=3:6
        push!(means,xf[1][j,(i-1)*3+2])
        push!(sems,xf[1][j,(i-1)*3+3])
    end    
    lines!(ax,t,means,color=getRandomColor(i))
    errorbars!(ax,t,means,sems,label=datalabel,color=getRandomColor(i))
end

legend = Legend(fig[1,2],ax)

display(fig)