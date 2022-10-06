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

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

fig = Figure(resolution=(6000,6000))

for (i,r) in enumerate(runs)
    maskData = load(datadir("exp_pro","masks",r[1:end-4],"$(r[1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData

    centroidData = load(datadir("exp_pro","emCentroidsInteractive","$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData

    # Import image file and convert to grayscale
    imageIn = load(fileName)
    imSize = size(imageIn)
    grayImage = Gray.(imageIn)

    ax2 = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect())
    image!(ax2,rotr90(imageIn))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=4)
    hidedecorations!(ax2)
    hidespines!(ax2)
end 


resize_to_layout!(fig)
display(fig)

save(datadir("exp_pro","emCentroidsInteractive","grid.png"),fig)



