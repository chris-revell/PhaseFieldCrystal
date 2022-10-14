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
using Statistics

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

mkpath(datadir("exp_pro","emCentroidMeasurements"))

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

fig = Figure(resolution=(6000,6000),backgroundcolor=:white,fontsize=64)

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
lengthPerPixelDict = Dict()
for r in runs 
    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
    lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
end

axes = Dict()
sizes = Dict()

for (i,r) in enumerate(runs)

    imageIn = load(datadir("exp_pro","cropped",r))
    
    centroidData = load(datadir("exp_pro","emCentroidsInteractive","$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(abs.([xs ys]))/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor

    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs,ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

    pairs = Tuple[]
    for t in eachrow(tri)
        for i=1:3
            labelSizeOrder = sort(t)
            (labelSizeOrder[1],labelSizeOrder[2]) in pairs ? nothing : push!(pairs,(labelSizeOrder[1],labelSizeOrder[2]))
            (labelSizeOrder[1],labelSizeOrder[3]) in pairs ? nothing : push!(pairs,(labelSizeOrder[1],labelSizeOrder[3]))
            (labelSizeOrder[2],labelSizeOrder[3]) in pairs ? nothing : push!(pairs,(labelSizeOrder[2],labelSizeOrder[3]))
        end
    end 

    lengths = norm.(centroidLocations[first.(pairs)]-centroidLocations[last.(pairs)])

    lengths .*= lengthPerPixelDict[r]*1000.0

    display(r)
    display(lengthPerPixelDict[r])
    display(mean(lengths))

    # lengths = Float64[]
    # for p in pairs 
    #     push!(lengths, norm(p[1]-p[2]))
    # end

    ax2 = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect(),backgroundcolor=:white)
    image!(ax2,rotr90(imageIn))
    #scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=ceil(Int64,10000/imSize[1]))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=10)
    for p in pairs
        lines!(ax2,Point.(centroidLocations[[p[1],p[2]]]))
    end
    hidedecorations!(ax2)
    hidespines!(ax2)
    Label(fig[(i-1)%6+1,(i-1)÷6+1, Bottom()], L"\mu=%$(round(mean(lengths),digits=1))nm,~\sigma=%$(round(std(lengths),digits=1))", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0))
    
    axes[r] = ax2
    sizes[r] = size(imageIn)
end 

lengthDict = Dict()
for r in runs 
    lengthDict[r] = lengthPerPixelDict[r].*sizes[r]
end

xMax = maximum(first.(values(lengthDict)))
yMax = maximum(last.(values(lengthDict)))

for r in runs 
    xlims!(axes[r],(0,xMax)./lengthPerPixelDict[r])
    ylims!(axes[r],(0,yMax)./lengthPerPixelDict[r])
end

resize_to_layout!(fig)
display(fig)

save(datadir("exp_pro","emCentroidMeasurements","grid.png"),fig)
