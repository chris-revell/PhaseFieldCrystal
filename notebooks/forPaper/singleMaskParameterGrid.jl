using DrWatson; @quickactivate
using Images
using CairoMakie
using Colors
using GeometryBasics
using LinearAlgebra
using FromFile
using FileIO
using ConcaveHull
using VoronoiCells
using GR: delaunay
using CSV
using DataFrames
using DelimitedFiles

@from "$(srcdir("ColourFunctions.jl"))" using ColourFunctions

# mask = "17tailT_4800X_HUI_0001_0"
# runs = [r for r in Vector(readdlm(datadir("exp_pro", "filesToUse.txt"))[:, 1]) if !(occursin("mp13ko", r) || occursin("18tailT_4800X_HUI_0007_0", r) || occursin("18tailT_4800X_HUI_0008_0", r))]

runs = [r for r in Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1]) if !(occursin("mp13ko",r) || occursin("18tailT_4800X_HUI_0007_0",r) || occursin("18tailT_4800X_HUI_0008_0",r) )]


for r in runs
    mask = r[1:end-4]

    # Collate results as a dataframe 
    results = collect_results(datadir("fromCSF","allMasksPhasespaceSeparateLengths",mask); subfolders = false)

    fig = Figure(resolution=(6000,6000),fontsize=64)

    axesDict = Dict()

    maskIn = load(datadir("exp_pro","masksCompressed",mask,"$mask.png"))
    maskImage = fill(RGBA(1,1,1,1),size(maskIn))
    for i=1:size(maskIn)[1]
        for j=1:size(maskIn)[2]
            if maskIn[i,j] < 0.5
                maskImage[i,j] = RGBA(0.0,0.0,0.0,1.0)
            else
                maskImage[i,j] = RGBA(0.0,0.0,0.0,0.0)
            end
        end
    end   

    # Loop to process data from each run 
    for i=1:nrow(results)
        display(i)
        # Convert simulation result to a 2D matrix
        uMat = reshape(results[i,:u][end],(results[i,:nY],results[i,:nX]))#transpose(reshape(results[i,:u],(results[i,:nY],results[i,:nX])))
        
        # # Concave hull to identify boundary fibrils 
        # hull = concave_hull(shiftedCentroidLocations,3)
        
        ax = CairoMakie.Axis(fig[mod(i-1,6)+1,(i-1)÷6+1],aspect=DataAspect())
        heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
        hidedecorations!(ax)
        hidespines!(ax)
        image!(ax,rotr90(maskImage))
        # Map parameters to axis in axes dictionary 
        axesDict[(results[i,:r],results[i,:ϕ0])] = ax
        
    end 

    # Rearrange axes to put relevant r and ϕ0 values in order
    sortedrs  = unique!(sort(first.(keys(axesDict))))
    sortedϕ0s = unique!(sort(last.(keys(axesDict))))
    for (i,r) in enumerate(sortedrs)
        for (j,ϕ0) in enumerate(sortedϕ0s)
            fig[length(sortedrs)+1-i,j] = axesDict[(r,ϕ0)]
            # Label(fig[length(sortedrs)+1-i,j,Bottom()],labelsDict[(r,ϕ0)])
        end
    end
    for (i,r) in enumerate(sortedrs)
        Label(fig[length(sortedrs)+1-i,1,Left()],"r=$(r)")
    end
    for (j,ϕ0) in enumerate(sortedϕ0s)
        Label(fig[end,j,Bottom()],"ϕ₀=$(ϕ0)")
    end

    # Resize columns
    for i=1:length(sortedϕ0s)
        colsize!(fig.layout,i,Aspect(1.0,1.0))
    end
    resize_to_layout!(fig)

    # display(fig)

    save(datadir("fromCSF","allMasksPhasespaceSeparateLengths",mask,"$(mask)ParameterGrid.png"),fig)
end


