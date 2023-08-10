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

@from "$(srcdir("ColourFunctions.jl"))" using ColourFunctions

runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])

fig = Figure(resolution=(6000,6000),fontsize=64,backgroundcolor=:white)

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

    ax2 = CairoMakie.Axis(fig[mod(i-1,6)+1,(i-1)÷6+1],aspect=DataAspect(),backgroundcolor=:white)
    image!(ax2,rotr90(imageIn))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=10)
    hidedecorations!(ax2)
    hidespines!(ax2)
    Label(fig[mod(i-1,6)+1,(i-1)÷6+1, BottomLeft()], "$i", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0))
    axes[r] = ax2
    sizes[r] = imSize
end 

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))

lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
lengthPerPixelDict = Dict()
lengthDict = Dict()
for r in runs 
    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
    lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
    lengthDict[r] = lengthPerPixelDict[r].*sizes[r]
end

xMax = maximum(first.(values(lengthDict)))
yMax = maximum(last.(values(lengthDict)))

for r in runs 
    xlims!(axes[r],(0,xMax)./lengthPerPixelDict[r])
    ylims!(axes[r],(0,yMax)./lengthPerPixelDict[r])
end

colgap!(fig.layout, 1, -700)
colgap!(fig.layout, 2, -600)
colgap!(fig.layout, 3, -600)
colgap!(fig.layout, 4, -500)
colgap!(fig.layout, 5, -100)

rowgap!(fig.layout, 1, -100)
rowgap!(fig.layout, 2, -100)
rowgap!(fig.layout, 3, -200)

resize_to_layout!(fig)
display(fig)

save(datadir("exp_pro","emCentroidsInteractive","emCentroidsGrid.png"),fig)



