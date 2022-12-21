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
@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions


function imageToMask(fileName,distanceGaussian,dilateCount,erodeCount,lX,h)
    
    mkpath(datadir("exp_pro","masks",splitpath(fileName)[end][1:end-4]))
    mkpath(datadir("exp_pro","masksCompressed",splitpath(fileName)[end][1:end-4]))
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
    
    # Create inter-cellular space mask
    # Prune segments below threshold size, adding pruned segments to the largest neighbouring segment
    size1 = 5000
    seg2 = prune_segments(seg1, i->(segment_pixel_count(seg1,i)<size1), (i,j)->(-segment_pixel_count(seg1,j)))
    seg3 = prune_segments(seg2, first.(sort(collect(seg2.segment_pixel_count), by=x->x[2])[1:end-2]), (i,j)->-segment_pixel_count(seg2,j))
    
    centralSegment = seg3.image_indexmap[imSize.÷2...]

    newIndexMap = copy(grayImage)
    for i=1:imSize[1]
        for j=1:imSize[2]
            if seg3.image_indexmap[i,j] == centralSegment
                newIndexMap[i,j] = Gray(1)
            else
                newIndexMap[i,j] = Gray(0)
            end
        end
    end

    safesave(datadir("exp_pro","masks",splitpath(fileName)[end][1:end-4],splitpath(fileName)[end]),newIndexMap)
    
    compressedNX = lX/h
    percentage_scale = compressedNX/size(imageIn,2)
    new_size = trunc.(Int, size(imageIn) .* percentage_scale)
    compressedMask = imresize(newIndexMap, new_size)
    for i in eachindex(compressedMask)
        if compressedMask[i]>0.5
            compressedMask[i]=Gray(1.0)
        else
            compressedMask[i]=Gray(0.0)
        end
    end
    safesave(datadir("exp_pro","masksCompressed",splitpath(fileName)[end][1:end-4],splitpath(fileName)[end]), compressedMask)

    fig = CairoMakie.Figure()
    ax1 = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax1)
    hidespines!(ax1)
    image!(ax1,rotr90(imageIn))
    ax2 = CairoMakie.Axis(fig[1,2],aspect=DataAspect())
    hidedecorations!(ax2)
    hidespines!(ax2)
    image!(ax2,rotr90(newIndexMap))
    safesave(datadir("exp_pro","masks",splitpath(fileName)[end][1:end-4],"$(splitpath(fileName)[end][1:end-4])_comparison.png"), fig)
    safesave(datadir("exp_pro","masks",splitpath(fileName)[end][1:end-4],"$(splitpath(fileName)[end][1:end-4]).jld2"),@strdict fileName distanceGaussian dilateCount erodeCount newIndexMap compressedMask lX h)

    return newIndexMap
end

runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])
distanceGaussian = 1.0
dilateCount = 2
erodeCount = 2
croppedLX = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","croppedLX.csv")))
lXDict = Dict(croppedLX[!,:file].=>croppedLX[!,:lX])
# okMasks   = [f for f in readdir(datadir("exp_pro","masks","ok","compressed")) if f[end-3:end]==".png"]
# subFrame  = filter(:file => f->f in okMasks, croppedLX)
scalingLX = 200.0/1.82314
for r in runs
    imageToMask(datadir("exp_pro","cropped",r),distanceGaussian,dilateCount,erodeCount,lXDict[r]*scalingLX,0.4)   
end
