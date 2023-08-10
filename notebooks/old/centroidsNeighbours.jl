using DrWatson; @quickactivate
using Images
using ImageBinarization
using ImageSegmentation
using Random
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

function neighbourColours(x)
    if x==6
        return :white
    elseif x==5 
        return :red 
    elseif x==7
        return :blue
    else 
        return :black 
    end
end

function centroidsToMeasurements(fileName)
    mkpath(datadir("exp_pro","emCentroidNeighbours"))
    
    maskData = load(datadir("exp_pro","masks",splitpath(fileName)[end][1:end-4],"$(splitpath(fileName)[end][1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData
    
    centroidData = load(datadir("exp_pro","emCentroidsInteractive","$(splitpath(fileName)[end][1:end-4]).jld2"))
    @unpack centroidLocations = centroidData
    
    # Import image file and convert to grayscale
    imageIn = load(fileName)
    imSize = size(imageIn)
    grayImage = Gray.(imageIn)

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum([xs ys])/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor

    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs,ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

    # Concave hull to identify boundary fibrils 
    hull = concave_hull(shiftedCentroidLocations,1)
    # Indices of fibrils within the hull 
    hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])

    # Voronoi tessellation of centroid positions within (0,0) (1,1) box
    tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))


    fig = Figure(resolution=(1000,500))
    ax1 = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    image!(ax1,rotr90(imageIn))
    scatter!(ax1,centroidLocations,color=(:orange,1.0),markersize=4)
    for i=1:n
        poly!(ax1,centroidLocations[tri[i,:]],color=(:white,0.0),strokecolor=(:orange,1.0),strokewidth=1.0)
    end

    ax2 = CairoMakie.Axis(fig[1,2],aspect=DataAspect())
    # Plot voronoi cells around each fibril coloured by neighbour count

    
    maskImage = fill(RGBA(1,1,1,1),size(imageIn))
    for i=1:size(imageIn)[1]
        for j=1:size(imageIn)[2]
            if newIndexMap[i,j] < 0.5
                maskImage[i,j] = RGBA(0.0,0.0,0.0,1.0)
            else
                maskImage[i,j] = RGBA(0.0,0.0,0.0,0.0)
            end
        end
    end    

 


    for (i,c) in enumerate(tess.Cells)
        if i ∉ hullInds
            poly!(ax2, c.*scalingFactor, color=neighbourColours(nNeighbours[i]),strokecolor=(:black,1.0),strokewidth=1.0)
        else
            # poly!(ax2, c.*scalingFactor, color=:white,strokecolor=(:black,1.0),strokewidth=1.0)
        end
    end
    image!(ax2,rotr90(maskImage))
    # poly!(ax2,hull.vertices,color=(:grey,1.0))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=4)
    hidedecorations!(ax2)
    hidespines!(ax2)
    resize_to_layout!(fig)
    display(fig)

    save(datadir("exp_pro","emCentroidNeighbours","$(splitpath(fileName)[end])"),fig)
end

# fileName = "/Users/christopher/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/cropped/mp13ko-3wiew_4800X_hui_0002_2.png"
runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])
# lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
# distanceGaussian = 2.0
# fibrilMinSize = 10
# dilateCount = 2
# erodeCount = 2
# fibrilMaxSize = 5000
# saveFlag=1
for r in runs
    centroidsToMeasurements(datadir("exp_pro","cropped",r))
end


