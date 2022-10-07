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

fig = Figure(resolution=(6000,6000),fontsize=64)

axes = Dict()
sizes = Dict()

for (i,r) in enumerate(runs)
    maskData = load(datadir("exp_pro","masks",r[1:end-4],"$(r[1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData

    centroidData = load(datadir("exp_pro","emCentroidsInteractive","$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData

    # Import image file and convert to grayscale
    imageIn = load(datadir("exp_pro","cropped",r))
    imSize = size(imageIn)
    grayImage = Gray.(imageIn)

    ax2 = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect())
    image!(ax2,rotr90(imageIn))
    #scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=ceil(Int64,10000/imSize[1]))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=10)
    hidedecorations!(ax2)
    hidespines!(ax2)
    Label(fig[(i-1)%6+1,(i-1)÷6+1, Bottom()], "$i", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0))
    axes[r] = ax2
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



