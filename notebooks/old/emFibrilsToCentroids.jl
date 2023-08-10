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
using DelimitedFiles
@from "$(srcdir("ColourFunctions.jl"))" using ColourFunctions


function emFibrilsToCentroids(fileName,distanceGaussian,dilateCount,erodeCount,fibrilMinSize,fibrilMaxSize,saveFlag)
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
    for i=1:erodeCount
        erode!(binarizedImage)
    end
    for i=1:dilateCount
        dilate!(binarizedImage)
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

    fig = Figure(resolution=(4000,4000))
    ax1 = CairoMakie.Axis(fig[2,1],aspect=DataAspect())
    ax2 = CairoMakie.Axis(fig[2,2],aspect=DataAspect())
    ax3 = CairoMakie.Axis(fig[1,2],aspect=DataAspect())
    ax4 = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax1); hidespines!(ax1)
    hidedecorations!(ax2); hidespines!(ax2)
    hidedecorations!(ax3); hidespines!(ax3)
    hidedecorations!(ax4); hidespines!(ax4)
    image!(ax1,rotr90(map(i->get_random_color(i), fibrilIndexMap)))
    image!(ax3,rotr90(binarizedImage))
    image!(ax4,rotr90(grayImage))
    image!(ax2,rotr90(imageIn))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=30)
    resize_to_layout!(fig)
    display(fig)
    
    if saveFlag==1
        safesave(datadir("exp_pro","emCentroids",splitpath(fileName)[end][1:end-4],splitpath(fileName)[end]),fig)
        safesave(datadir("exp_pro","emCentroids",splitpath(fileName)[end][1:end-4],"$(splitpath(fileName)[end][1:end-4]).jld2"),@strdict fileName fibrilMinSize fibrilMaxSize distanceGaussian dilateCount erodeCount centroidLocations fig)
    end 

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

