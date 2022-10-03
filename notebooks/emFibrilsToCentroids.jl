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
using CSV
using DataFrames
@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions


function emFibrilsToCentroids(fileName,distanceGaussian,dilateCount,erodeCount,fibrilMinSize,fibrilMaxSize)
    mkpath(datadir("exp_pro","emCentroids",splitpath(fileName)[end][1:end-4]))
    
    # Import previously generated mask
    maskData = load(datadir("exp_pro","masks",splitpath(fileName)[end][1:end-4],"$(splitpath(fileName)[end][1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData
    
    # Import image file and convert to grayscale
    imageIn = load(fileName)
    imSize = size(imageIn)
    grayImage = Gray.(imageIn)
    # Apply Gaussian filter
    filteredImage  = imfilter(grayImage,Kernel.gaussian(distanceGaussian))
    # Binarise with Otsu algorithm
    binarizedImage = binarize(filteredImage,Otsu())
    # Erode and dilate    
    for i=1:dilateCount
        dilate!(binarizedImage)
    end
    for i=1:erodeCount
        erode!(binarizedImage)
    end

    # Initial segmentation 
    seg1 = fast_scanning(binarizedImage, 0.01)
    
    # Exclude anything not in the intra-cellular space from the segmentation using the imported mask
    centralSegment = newIndexMap[imSize.÷2...]
    fibrilIndexMap = copy(seg1.image_indexmap)
    for i=1:imSize[1]
        for j=1:imSize[2]
            if newIndexMap[i,j] != centralSegment
                fibrilIndexMap[i,j] = 0
            end
        end
    end

    # Find segment labels within fibril index map
    labels = unique(fibrilIndexMap)
    # Remove 0 segment, used to represent intra-cell space
    filter!(val->val≠0,labels)
    centroidLocations = Point2[]
    for k in labels
        pixels = zeros(2)
        count = 0
        for i=1:imSize[1]
            for j=1:imSize[2]
                if fibrilIndexMap[i,j] == k
                    pixels .+= [j,-i]
                    count += 1
                end
            end
        end
        # Filter objects smaller than fibrilMinSize and larger than some threshold big enough to only exclude the background segments. 
        fibrilMinSize<count<fibrilMaxSize ? push!(centroidLocations,Point2(pixels./count)) : nothing
    end
    centroidLocations .+= Point2(0.0,imSize[1]*1.0)

    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    image!(ax,rotr90(imageIn))
    scatter!(ax,centroidLocations,color=(:orange,1.0),markersize=6)
    safesave(datadir("exp_pro","emCentroids",splitpath(fileName)[end][1:end-4],splitpath(fileName)[end]),fig)
    safesave(datadir("exp_pro","emCentroids",splitpath(fileName)[end][1:end-4],"$(splitpath(fileName)[end][1:end-4]).jld2"),@strdict fileName fibrilMinSize fibrilMaxSize distanceGaussian dilateCount erodeCount centroidLocations fig)
    # display(fig)
    return fig, centroidLocations, fibrilIndexMap
end

# fileName = "/Users/christopher/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/cropped/mp13ko-3wiew_4800X_hui_0002_2.png"
runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])
# lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
distanceGaussian = 1.0
fibrilMinSize = 10
fibrilMaxSize = 1000
dilateCount = 2
erodeCount = 2
for r in runs
    emFibrilsToCentroids(datadir("exp_pro","cropped",r),distanceGaussian,dilateCount,erodeCount,fibrilMinSize,fibrilMaxSize)
end

