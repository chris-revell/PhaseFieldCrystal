using DrWatson
@quickactivate
using UnPack
using CairoMakie; set_theme!(figure_padding=0,fontsize=32)
using ColorSchemes
using DifferentialEquations
using JLD2
using LaTeXStrings
using DataFrames
using Images
using ConcaveHull

# Location of data within data/fromCSF/ 
folderPath = "PhaseSpace4"

# Collate results as a dataframe 
results = collect_results(datadir("fromCSF",folderPath); subfolders = true)

# Set up figure canvas and dictionary to map parameters to axes
fig = Figure(resolution=(6000,6000))
axesDict = Dict()

# Loop to process data from each run 
for i=1:nrow(subset(results, :m => m -> m.== 0.1))
    display(i)
    # Convert simulation result to a 2D matrix
    uMat = reshape(results[i,:u],(results[i,:nY],results[i,:nX]))#transpose(reshape(results[i,:u],(results[i,:nY],results[i,:nX])))
    
    # # Concave hull to identify boundary fibrils 
    # hull = concave_hull(shiftedCentroidLocations,3)
    
    ax = CairoMakie.Axis(fig[i%6+1,i÷6+1],aspect=DataAspect())
    heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
    hidedecorations!(ax)
    hidespines!(ax)
    # Map parameters to axis in axes dictionary 
    axesDict[(results[i,:r],results[i,:ϕ0])] = ax
    
end 

# Rearrange axes to put relevant r and ϕ0 values in order
sortedrs  = unique!(sort(first.(keys(axes))))
sortedϕ0s = unique!(sort(last.(keys(axes))))
for (i,r) in enumerate(sortedrs)
    for (j,ϕ0) in enumerate(sortedϕ0s)
        fig[length(sortedrs)+1-i,length(sortedϕ0s)+1-j][1,1] = axesDict[(r,ϕ0)]
        # Colorbar(fig[length(sortedrs)-i,length(sortedϕ0s)-j][1,2], limits=(-1,1),colormap=:bwr)#, size = 25)
        Label(fig[length(sortedrs)+1-i,length(sortedϕ0s)+1-j,Bottom()],L"r=%$r, \phi_0=%$ϕ0",fontsize=64)
    end
end

# Resize columns
for i=1:length(sortedϕ0s)
    colsize!(fig.layout,i,Aspect(1.0,1.0))
end
resize_to_layout!(fig)

display(fig)

save(datadir("fromCSF",folderPath,"resultPhaseSpace.png"),fig)
