using XLSX
using DataFrames
using StatsBase
using CairoMakie
using FromFile

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

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

fullGeneNames = Vector(TailTendon_TimePoints_MSqRob2_JC[1]["B2:B2887"][:,1])

indices = findall(x->(!ismissing(x)&&x∈subsetGeneNames), fullGeneNames)
counts = countmap(fullGeneNames[indices])
# for key in keys(counts)
#     if counts[key]>1
#         display(key)
#     end
# end
indicesToRemove = Int64[]
for i in 1:length(indices)
    if counts[fullGeneNames[indices[i]]]>1 && ismissing(TailTendon_TimePoints_MSqRob2_JC[1][indices[i]+1,3])
        push!(indicesToRemove,i)
        # display(fullGeneNames[indices[i]])
        # display(TailTendon_TimePoints_MSqRob2_JC[1][indices[i]+1,3])
    end
end
deleteat!(indices,indicesToRemove)

for (sheet,day) in enumerate(names(foldChangeDataFrame)[4:end])
    data = TailTendon_TimePoints_MSqRob2_JC[sheet]["C2:D2887"]
    foldChangeDataFrame[!,day] .= data[indices,1]
    semDataFrame[!,day] .= data[indices,2]
end

XLSX.openxlsx(datadir("exp_pro","MSSpreadsheets","foldChange.xlsx"),mode="w") do xf
    # rename!(xf[1], "foldChanges")
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
    save(datadir("exp_pro","MSSpreadsheets","$l fold change.png"),fig)
end

for l in differentLabels
    fig = Figure(resolution=(500,500))    
    subsetFoldChange = dropmissing(filter(:label=>n->n==l,foldChangeDataFrame))
    proteins = Vector(subsetFoldChange[!,2])
    heatmapData = Matrix(subsetFoldChange[!,3:7])
    ax = Axis(fig[1,1], xticks=(1:5, string.(t)), yticks=(1:length(proteins), reverse(proteins)))#, aspect=DataAspect())
    heatmap!(ax,rotr90(heatmapData), colormap=:bwr, colorrange=(-maximum(heatmapData)+1, maximum(heatmapData)+1))    
    ax.xlabel = "Time /days"
    ax.ylabel = "Protein"
    ax.title = l  
    Colorbar(fig[1,2], limits=(-maximum(heatmapData)+1, maximum(heatmapData)+1), colormap=:bwr)  
    # resize_to_layout!(fig)
    save(datadir("exp_pro","MSSpreadsheets","$l fold change heatmap.png"),fig)
end

fig = Figure(resolution=(500,2000))
allFoldChange = dropmissing(foldChangeDataFrame)
allProteins = Vector(allFoldChange[!,2])
allHeatmapData = Matrix(allFoldChange[!,3:7])
ax = Axis(fig[1:5,1], xticks=(1:5, string.(t)), yticks=(1:length(allProteins), reverse(allProteins)))#, aspect=DataAspect())
heatmap!(ax,rotr90(allHeatmapData), colormap=:bwr, colorrange=(-maximum(allHeatmapData)+1, maximum(allHeatmapData)+1))    
ax.xlabel = "Time /days"
ax.ylabel = "Protein"
ax.title = "All proteins"
Colorbar(fig[2:4,2], limits=(-maximum(allHeatmapData)+1, maximum(allHeatmapData)+1), colormap=:bwr, label="Fold change")  
resize_to_layout!(fig)
save(datadir("exp_pro","MSSpreadsheets","All proteins fold change heatmap.png"),fig)