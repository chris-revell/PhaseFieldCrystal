# This script assumes relevant spreadsheets are found in data/exp_pro/MSSpreadsheets
using Statistics
using XLSX
using DataFrames
using StatsBase
using CairoMakie
using FromFile

@from "$(srcdir("ColourFunctions.jl"))" using ColourFunctions

Values_for_heatmap = XLSX.readxlsx(datadir("exp_pro","MSSpreadsheets","Values for heatmap.xlsx"))[1]

subsetGeneNames = Values_for_heatmap["B2:B107"][:,1]
geneLabels = Values_for_heatmap["A2:A107"][:,1]

TailTendon_TimePoints_MSqRob2_JC = XLSX.readxlsx(datadir("exp_pro","MSSpreadsheets","TailTendon_TimePoints_MSqRob2_JC.xlsx"))

foldChangeDataFrame = DataFrame(
    label = geneLabels,
    gene_name = subsetGeneNames, 
    E12_5 = ones(Float64,length(subsetGeneNames)),
    E13_0 = zeros(Float64,length(subsetGeneNames)),
    E13_5 = zeros(Float64,length(subsetGeneNames)),
    E14_0 = zeros(Float64,length(subsetGeneNames)),
    E14_5 = zeros(Float64,length(subsetGeneNames))
)

semDataFrame = DataFrame(
    label = geneLabels,
    gene_name = subsetGeneNames, 
    E12_5 = ones(Float64,length(subsetGeneNames)),
    E13_0 = zeros(Float64,length(subsetGeneNames)),
    E13_5 = zeros(Float64,length(subsetGeneNames)),
    E14_0 = zeros(Float64,length(subsetGeneNames)),
    E14_5 = zeros(Float64,length(subsetGeneNames))
)

fullGeneNames = Vector(TailTendon_TimePoints_MSqRob2_JC[1]["B1:B2887"][:,1])

indices = Int64[]
for gene in subsetGeneNames
    index = findall(x->(!ismissing(x)&&x==gene), fullGeneNames)
    append!(indices,index)
end
counts = countmap(fullGeneNames[indices])

indicesToRemove = Int64[]
for i in 1:length(indices)
    if counts[fullGeneNames[indices[i]]]>1 && ismissing(TailTendon_TimePoints_MSqRob2_JC[1][indices[i]+1,3])
        push!(indicesToRemove,i)
    end
end
deleteat!(indices,indicesToRemove)

for (sheet,day) in enumerate(names(foldChangeDataFrame)[4:end])
    data = TailTendon_TimePoints_MSqRob2_JC[sheet]["C1:D2887"]
    foldChangeDataFrame[!,day] .= data[indices,1]
    semDataFrame[!,day] .= data[indices,2]
end

XLSX.openxlsx(datadir("exp_pro","MSSpreadsheets","foldChange.xlsx"),mode="w") do xf
    XLSX.addsheet!(xf,"foldChanges")
    sheet = xf["foldChanges"]
    for i in 1:length(names(foldChangeDataFrame))
        sheet[XLSX.CellRef(1, i)] = names(foldChangeDataFrame)[i]
    end
    for r in 1:size(foldChangeDataFrame,1), c in 1:size(foldChangeDataFrame,2)
         sheet[XLSX.CellRef(r+1, c)] = foldChangeDataFrame[r,c]
    end
    XLSX.addsheet!(xf,"sem")
    sheet = xf["sem"]
    for i in 1:length(names(foldChangeDataFrame))
        sheet[XLSX.CellRef(1, i)] = names(foldChangeDataFrame)[i]
    end
    for r in 1:size(semDataFrame,1), c in 1:size(semDataFrame,2)
         sheet[XLSX.CellRef(r+1, c)] = semDataFrame[r,c]
    end
end

differentLabels = unique(geneLabels)
t = [12.5,13.0,13.5,14.0,14.5]
for l in differentLabels
    fig = Figure(resolution=(500,500))
    ax = Axis(fig[1,1], xgridvisible = false, ygridvisible = false)
    subsetFoldChange = filter(:label=>n->n==l,foldChangeDataFrame)
    subsetSem = filter(:label=>n->n==l,semDataFrame)
    for i=1:nrow(subsetFoldChange)
        if true ∉ ismissing.(Vector(subsetFoldChange[i,3:end]))# || true ∉ ismissing.(Vector(subsetSem[i,3:end]))
            lines!(ax,t,Vector(subsetFoldChange[i,3:end]),label=subsetFoldChange[i,2],color=get_random_color(i), linewidth=2)
            errorbars!(ax,t,Vector(subsetFoldChange[i,3:end]),[0.0, Vector(subsetSem[i,4:end])...],color=get_random_color(i), linewidth=2)
        end
    end
    legend = Legend(fig[1,2],ax)
    ax.xlabel = "Time /days"
    ax.ylabel = "Fold change"
    ax.title = l
    xlims!(12.5, 14.6)
    resize_to_layout!(fig)
    #save(datadir("exp_pro","MSSpreadsheets","$l fold change.png"),fig)
end

for l in differentLabels
    fig = Figure(resolution=(500,500))    
    subsetFoldChange = dropmissing(filter(:label=>n->n==l,foldChangeDataFrame))
    proteins = Vector(subsetFoldChange[!,2])
    proteins[proteins.=="P3h1;Lepre1"] .= "P3h1"
    heatmapData = Matrix(subsetFoldChange[!,3:7])
    ax = Axis(fig[1,1], xticks=(1:5, string.(t)), yticks=(1:length(proteins), reverse(proteins)))#, aspect=DataAspect())
    heatmap!(ax,rotr90(heatmapData), colormap=:bwr, colorrange=(-maximum(abs.(heatmapData))+1, maximum(abs.(heatmapData))+1))    
    ax.xlabel = "Time /days"
    ax.ylabel = "Protein"
    ax.title = l  
    Colorbar(fig[1,2], limits=(-maximum(heatmapData)+1, maximum(heatmapData)+1), colormap=:bwr, label="Fold change")  
    # resize_to_layout!(fig)
    #save(datadir("exp_pro","MSSpreadsheets","$l fold change heatmap.png"),fig)
end

fig = Figure(resolution=(1500,1000))
allFoldChange = dropmissing(foldChangeDataFrame)
allProteins = Vector(allFoldChange[!,2])
allProteins[allProteins.=="P3h1;Lepre1"] .= "P3h1"
allProteinsRearranged = copy(allProteins)
allProteinsRearranged[1:25] .= reverse(allProteins[1:25])
allProteinsRearranged[26:50] .= reverse(allProteins[26:50])
allProteinsRearranged[51:75] .= reverse(allProteins[51:75])
allHeatmapData = Matrix(allFoldChange[!,3:7])
allLabels = Vector(allFoldChange[!,1])
numericalLabels = [findall(x->x==f,unique(allLabels))[1] for f in allLabels]
colors = [:black, :red, :orange, :green, :blue, :indigo, :violet]
colorLabels = colors[numericalLabels]
colorLabelsRearranged = copy(colorLabels)
colorLabelsRearranged[1:25] .= reverse(colorLabels[1:25])
colorLabelsRearranged[26:50] .= reverse(colorLabels[26:50])
colorLabelsRearranged[51:75] .= reverse(colorLabels[51:75])

ax1 = Axis(fig[1,1], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(1:25, colorLabelsRearranged[1:25])]))
ax2 = Axis(fig[1,2], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(26:50, colorLabelsRearranged[26:50])]))
ax3 = Axis(fig[1,3], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(51:75, colorLabelsRearranged[51:75])]))
heatmap!(ax1,rotr90(allHeatmapData[1:25,:]), colormap=:bwr, colorrange=(-maximum(abs.(allHeatmapData))+1, maximum(abs.(allHeatmapData))+1))    
heatmap!(ax2,rotr90(allHeatmapData[26:50,:]), colormap=:bwr, colorrange=(-maximum(abs.(allHeatmapData))+1, maximum(abs.(allHeatmapData))+1))
heatmap!(ax3,rotr90(allHeatmapData[51:75,:]), colormap=:bwr, colorrange=(-maximum(abs.(allHeatmapData))+1, maximum(abs.(allHeatmapData))+1))
grid = GridLayout(fig[2,:])
# Box(grid[1,:], color = :black)
# Box(gd[i, 3], color = :gray90)
Label(grid[1, 1], justification = :center, "Protein colour labels:", fontsize=16, lineheight = 2.0)
Label(grid[1, 1+1], justification = :center, "$(unique(allLabels)[1])", color=colors[1], fontsize=16, lineheight = 2.0)
Label(grid[1, 2+1], justification = :center, "$(unique(allLabels)[2])", color=colors[2], fontsize=16, lineheight = 2.0)
Label(grid[1, 3+1], justification = :center, "$(unique(allLabels)[3])", color=colors[3], fontsize=16, lineheight = 2.0)
Label(grid[1, 4+1], justification = :center, "$(unique(allLabels)[4])", color=colors[4], fontsize=16, lineheight = 2.0)
Label(grid[1, 5+1], justification = :center, "$(unique(allLabels)[5])", color=colors[5], fontsize=16, lineheight = 2.0)
Label(grid[1, 6+1], justification = :center, "$(unique(allLabels)[6])", color=colors[6], fontsize=16, lineheight = 2.0)
Label(grid[1, 7+1], justification = :center, "$(unique(allLabels)[7])", color=colors[7], fontsize=16, lineheight = 2.0)
Box(grid[1,1:8], color=(:white,0.0))
ax2.xlabel = "Time /days"
ax1.ylabel = "Protein"
Colorbar(fig[1,4], limits=(-maximum(allHeatmapData)+1, maximum(allHeatmapData)+1), colormap=:bwr, label="Fold change") 

rowsize!(fig.layout,2,0.2)
resize_to_layout!(fig)
display(fig)
#save(datadir("exp_pro","MSSpreadsheets","All proteins fold change heatmap2.png"),fig)
#save(datadir("exp_pro","MSSpreadsheets","All proteins fold change heatmap2.svg"),fig)




allProteinsReference = Vector(allFoldChange[!,2])

absoluteChangeDataFrameLogsOriginal = DataFrame(
    label = geneLabels,
    gene_name = subsetGeneNames, 
    E12_5 = zeros(Float64,length(subsetGeneNames)),
    E13_0 = zeros(Float64,length(subsetGeneNames)),
    E13_5 = zeros(Float64,length(subsetGeneNames)),
    E14_0 = zeros(Float64,length(subsetGeneNames)),
    E14_5 = zeros(Float64,length(subsetGeneNames))
)

absoluteChangeDataFrameLogsOriginal[!,3:7] = Values_for_heatmap["C2:G107"]

absoluteChangeDataFrameOriginal = DataFrame(
    label = geneLabels,
    gene_name = subsetGeneNames, 
    E12_5 = zeros(Float64,length(subsetGeneNames)),
    E13_0 = zeros(Float64,length(subsetGeneNames)),
    E13_5 = zeros(Float64,length(subsetGeneNames)),
    E14_0 = zeros(Float64,length(subsetGeneNames)),
    E14_5 = zeros(Float64,length(subsetGeneNames))
)

allProteinsReference = copy(allProteins)

absoluteChangeDataFrameOriginal[!,3:7] = 2.0.^Values_for_heatmap["C2:G107"]

absoluteChangeDataFrame = filter(:gene_name=>n->n∈allProteinsReference, absoluteChangeDataFrameOriginal)
absoluteChangeDataFrameLogs = filter(:gene_name=>n->n∈allProteinsReference, absoluteChangeDataFrameLogsOriginal)

fig2 = Figure(resolution=(1500,1000))
allAbsoluteChange = dropmissing(absoluteChangeDataFrame)
allProteins = Vector(allAbsoluteChange[!,2])
allProteins[allProteins.=="P3h1;Lepre1"] .= "P3h1"
allProteinsRearranged = copy(allProteins)
allProteinsRearranged[1:25] .= reverse(allProteins[1:25])
allProteinsRearranged[26:50] .= reverse(allProteins[26:50])
allProteinsRearranged[51:75] .= reverse(allProteins[51:75])
allHeatmapData = Matrix(allAbsoluteChange[!,3:7])
allLabels = Vector(allAbsoluteChange[!,1])
numericalLabels = [findall(x->x==f,unique(allLabels))[1] for f in allLabels]
colors = [:black, :red, :orange, :green, :blue, :indigo, :violet]
colorLabels = colors[numericalLabels]
colorLabelsRearranged = copy(colorLabels)
colorLabelsRearranged[1:25] .= reverse(colorLabels[1:25])
colorLabelsRearranged[26:50] .= reverse(colorLabels[26:50])
colorLabelsRearranged[51:75] .= reverse(colorLabels[51:75])

ax1 = Axis(fig2[1,1], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(1:25, colorLabelsRearranged[1:25])]))
ax2 = Axis(fig2[1,2], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(26:50, colorLabelsRearranged[26:50])]))
ax3 = Axis(fig2[1,3], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(51:75, colorLabelsRearranged[51:75])]))
heatmap!(ax1,rotr90(allHeatmapData[1:25,:]), colormap=:batlow, colorrange=(minimum(allHeatmapData), maximum(allHeatmapData)))    
heatmap!(ax2,rotr90(allHeatmapData[26:50,:]), colormap=:batlow, colorrange=(minimum(allHeatmapData), maximum(allHeatmapData)))
heatmap!(ax3,rotr90(allHeatmapData[51:75,:]), colormap=:batlow, colorrange=(minimum(allHeatmapData), maximum(allHeatmapData)))
grid = GridLayout(fig2[2,:])
# Box(grid[1,:], color = :black)
# Box(gd[i, 3], color = :gray90)
Label(grid[1, 1], justification = :center, "Protein colour labels:", fontsize=16, lineheight = 2.0)
Label(grid[1, 1+1], justification = :center, "$(unique(allLabels)[1])", color=colors[1], fontsize=16, lineheight = 2.0)
Label(grid[1, 2+1], justification = :center, "$(unique(allLabels)[2])", color=colors[2], fontsize=16, lineheight = 2.0)
Label(grid[1, 3+1], justification = :center, "$(unique(allLabels)[3])", color=colors[3], fontsize=16, lineheight = 2.0)
Label(grid[1, 4+1], justification = :center, "$(unique(allLabels)[4])", color=colors[4], fontsize=16, lineheight = 2.0)
Label(grid[1, 5+1], justification = :center, "$(unique(allLabels)[5])", color=colors[5], fontsize=16, lineheight = 2.0)
Label(grid[1, 6+1], justification = :center, "$(unique(allLabels)[6])", color=colors[6], fontsize=16, lineheight = 2.0)
Label(grid[1, 7+1], justification = :center, "$(unique(allLabels)[7])", color=colors[7], fontsize=16, lineheight = 2.0)
Box(grid[1,1:8], color=(:white,0.0))
ax2.xlabel = "Time /days"
ax1.ylabel = "Protein"
Colorbar(fig2[1,4], limits=(minimum(allHeatmapData), maximum(allHeatmapData)), colormap=:batlow, label="Protein intensity") 

rowsize!(fig2.layout,2,0.2)
resize_to_layout!(fig2)
display(fig2)
save(datadir("exp_pro","MSSpreadsheets","All proteins absolute change heatmap.png"),fig2)
#save(datadir("exp_pro","MSSpreadsheets","All proteins fold change heatmap2.svg"),fig)








fig3 = Figure(resolution=(1500,1000))
allAbsoluteChange = dropmissing(absoluteChangeDataFrameLogs)
allProteins = Vector(allAbsoluteChange[!,2])
allProteins[allProteins.=="P3h1;Lepre1"] .= "P3h1"
allProteinsRearranged = copy(allProteins)
allProteinsRearranged[1:25] .= reverse(allProteins[1:25])
allProteinsRearranged[26:50] .= reverse(allProteins[26:50])
allProteinsRearranged[51:75] .= reverse(allProteins[51:75])
allHeatmapData = Matrix(allAbsoluteChange[!,3:7])
allLabels = Vector(allAbsoluteChange[!,1])
numericalLabels = [findall(x->x==f,unique(allLabels))[1] for f in allLabels]
colors = [:black, :red, :orange, :green, :blue, :indigo, :violet]
colorLabels = colors[numericalLabels]
colorLabelsRearranged = copy(colorLabels)
colorLabelsRearranged[1:25] .= reverse(colorLabels[1:25])
colorLabelsRearranged[26:50] .= reverse(colorLabels[26:50])
colorLabelsRearranged[51:75] .= reverse(colorLabels[51:75])

ax1 = Axis(fig3[1,1], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(1:25, colorLabelsRearranged[1:25])]))
ax2 = Axis(fig3[1,2], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(26:50, colorLabelsRearranged[26:50])]))
ax3 = Axis(fig3[1,3], xticks=(1:5, string.(t)), yticks = (1:25, [rich(allProteinsRearranged[val]; color) for (val, color) in zip(51:75, colorLabelsRearranged[51:75])]))
heatmap!(ax1,rotr90(allHeatmapData[1:25,:]), colormap=:batlow, colorrange=(minimum(allHeatmapData), maximum(allHeatmapData)))    
heatmap!(ax2,rotr90(allHeatmapData[26:50,:]), colormap=:batlow, colorrange=(minimum(allHeatmapData), maximum(allHeatmapData)))
heatmap!(ax3,rotr90(allHeatmapData[51:75,:]), colormap=:batlow, colorrange=(minimum(allHeatmapData), maximum(allHeatmapData)))
grid = GridLayout(fig3[2,:])
# Box(grid[1,:], color = :black)
# Box(gd[i, 3], color = :gray90)
Label(grid[1, 1], justification = :center, "Protein colour labels:", fontsize=16, lineheight = 2.0)
Label(grid[1, 1+1], justification = :center, "$(unique(allLabels)[1])", color=colors[1], fontsize=16, lineheight = 2.0)
Label(grid[1, 2+1], justification = :center, "$(unique(allLabels)[2])", color=colors[2], fontsize=16, lineheight = 2.0)
Label(grid[1, 3+1], justification = :center, "$(unique(allLabels)[3])", color=colors[3], fontsize=16, lineheight = 2.0)
Label(grid[1, 4+1], justification = :center, "$(unique(allLabels)[4])", color=colors[4], fontsize=16, lineheight = 2.0)
Label(grid[1, 5+1], justification = :center, "$(unique(allLabels)[5])", color=colors[5], fontsize=16, lineheight = 2.0)
Label(grid[1, 6+1], justification = :center, "$(unique(allLabels)[6])", color=colors[6], fontsize=16, lineheight = 2.0)
Label(grid[1, 7+1], justification = :center, "$(unique(allLabels)[7])", color=colors[7], fontsize=16, lineheight = 2.0)
Box(grid[1,1:8], color=(:white,0.0))
ax2.xlabel = "Time /days"
ax1.ylabel = "Protein"
Colorbar(fig3[1,4], limits=(minimum(allHeatmapData), maximum(allHeatmapData)), colormap=:batlow, label="log₂(Protein intensity)") 

rowsize!(fig3.layout,2,0.2)
resize_to_layout!(fig3)
display(fig3)
save(datadir("exp_pro","MSSpreadsheets","logAbundancesHeatmap.png"),fig3)
#save(datadir("exp_pro","MSSpreadsheets","All proteins fold change heatmap2.svg"),fig)

allAbsoluteChange = dropmissing(absoluteChangeDataFrame)
genenames = String[]
means = Float64[]
for i in eachrow(allAbsoluteChange)
    push!(genenames,i[2])
    push!(means,mean(i[3:7]))
end
genenames[genenames.=="P3h1;Lepre1"] .= "P3h1"
numericalLabels = [findall(x->x==f,unique(Vector(allAbsoluteChange[!,1])))[1] for f in allLabels]
colors = [:black, :red, :orange, :green, :blue, :indigo, :violet]
colorLabels = colors[numericalLabels]

fig4 = Figure(resolution=(1500,1000),fontsize=20)
ax9 = Axis(fig4[1,1], xticks = (1:75, [rich(genenames[val]; color) for (val, color) in zip(1:75, colorLabels)]),xticklabelrotation=π/2)
barplot!(ax9,collect(1:length(genenames)),log2.(means),color=colorLabels)
ax9.xlabelsize = 32
ax9.xlabel = "Gene name"
ax9.ylabelsize = 32
ax9.ylabel = "log₂(Time averaged protein intensity)"
ylims!(ax9,(20,33))
xlims!(0.5,75.5)
display(fig4)
save(datadir("exp_pro","MSSpreadsheets","timeAveragedAbundance.png"),fig4)