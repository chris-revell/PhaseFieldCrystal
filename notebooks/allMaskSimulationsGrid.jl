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

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

# Collate results as a dataframe 
results = collect_results!(datadir("fromCSF","allMasks"); subfolders = true)

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

fig = Figure(resolution=(6000,6000),fontsize=64)

axes = Dict()
sizes = Dict()

for (i,r) in enumerate(runs)

    subsetResults = subset(results, :path => m -> occursin.(r[1:end-4],m))

    maskData = load(datadir("exp_pro","masks",r[1:end-4],"$(r[1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData
    maskImage = fill(RGBA(1,1,1,1),size(newIndexMap))
    for i=1:size(newIndexMap)[1]
        for j=1:size(newIndexMap)[2]
            if newIndexMap[i,j] < 0.5
                maskImage[i,j] = RGBA(0.0,0.0,0.0,1.0)
            else
                maskImage[i,j] = RGBA(0.0,0.0,0.0,0.0)
            end
        end
    end    

    ax = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect())
    uMat = reshape(subsetResults[1,:u],(subsetResults[1,:nY],subsetResults[1,:nX]))
    heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
    image!(ax2,rotr90(maskImage))
    hidedecorations!(ax)
    hidespines!(ax)
    Label(fig[(i-1)%6+1,(i-1)÷6+1, Bottom()], "$i", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0))
    axes[r] = ax
    sizes[r] = imSize
end 

xMax = maximum(last.(values(sizes)))
yMax = maximum(first.(values(sizes)))

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))

lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
lengthPerPixelDict = Dict()
for r in runs 
    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
    lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
end


for r in runs 
    xlims!(axes[r],((last(sizes[r])-xMax)/2,(xMax+last(sizes[r]))/2).*lengthPerPixelDict[r])
    ylims!(axes[r],((first(sizes[r])-yMax)/2,(yMax+first(sizes[r]))/2).*lengthPerPixelDict[r])
end

resize_to_layout!(fig)
display(fig)

save(datadir("exp_pro","emCentroidsInteractive","grid.png"),fig)



