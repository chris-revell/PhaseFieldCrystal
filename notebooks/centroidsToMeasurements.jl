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

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions


function centroidsToMeasurements(fileName,distanceGaussian,dilateCount,erodeCount,size1,saveFlag)
    mkpath(datadir("exp_pro","emCentroidMeasurements"))
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
    seg2 = prune_segments(seg1, i->(segment_pixel_count(seg1,i)<size1), (i,j)->(-segment_pixel_count(seg1,j)))
    seg3 = prune_segments(seg2, first.(sort(collect(seg2.segment_pixel_count), by=x->x[2])[1:end-2]), (i,j)->-segment_pixel_count(seg2,j))
    intraCellSpace = (sort(collect(seg3.segment_pixel_count), by=x->x[2]))[1]

    newIndexMap = copy(seg1.image_indexmap)
    for i=1:size(imageIn)[1]
        for j=1:size(imageIn)[2]
            if seg3.image_indexmap[i,j] == intraCellSpace.first
                newIndexMap[i,j] = 0
            end
        end
    end

    # centroidLocations = findCentroidLocations(newIndexMap,size(imageIn),extraCellSpace)
    # Find segment labels within index map
    labels = unique(newIndexMap)
    # Remove 0 segment, used to represent intra-cell space
    filter!(val->val≠0,labels)
    # Remove extraCellSpace segment, used to represent extra-cellular space
    # filter!(val->val≠extraCellSpace,labels)
    centroidLocations = Point2[]
    for k in labels
        pixels = zeros(2)
        count = 0
        for i=1:imSize[1]
            for j=1:imSize[2]
                if newIndexMap[i,j] == k
                    pixels .+= [j,-i]
                    count += 1
                end
            end
        end
        # Filter objects smaller than fibrilMinSize and larger than some threshold big enough to only exclude the background segments. 
        fibrilMinSize<count<10000 ? push!(centroidLocations,Point2(pixels./count)) : nothing
    end
    centroidLocations .+= Point2(0.0,imSize[1]*1.0)

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


    fig = Figure()
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
    poly!(ax2,hull.vertices,color=(:black,0.5))
    # Map neighbour counts to colours 
    cmap = cgrad(:bwr, 5, categorical=true, rev=true)
    for (i,c) in enumerate(tess.Cells)
        # if i ∉ hullInds
            poly!(ax2, c, color=[nNeighbours[i]], colormap=cmap, colorrange=(4,8),highclip=:black,lowclip=:black)
        # end
    end
    hidedecorations!(ax2)
    hidespines!(ax2)
    xlims!(ax2,0,size(imageIn)[2]/scalingFactor)
    ylims!(ax2,0,size(imageIn)[1]/scalingFactor)
    

    if saveFlag == 1
        save(datadir("exp_pro","emCentroids","$(splitpath(fileName)[end])"),fig)
        save(datadir("exp_pro","emCentroids","$(splitpath(fileName)[end][1:end-4]).jld2"),@strdict fileName fibrilMinSize distanceGaussian dilateCount erodeCount centroidLocations fig)
    end
    display(fig)
    display(maximum(nNeighbours))
    display(minimum(nNeighbours))
    
    return fig, centroidLocations, newIndexMap
end

# fileName = "/Users/christopher/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/cropped/mp13ko-3wiew_4800X_hui_0002_2.png"
runs = [f for f in readdir(datadir("exp_pro","cropped")) if f[end-3:end]==".png"]
lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
distanceGaussian = 1.0
fibrilMinSize = 250
dilateCount = 2
erodeCount = 2
size1 = 5000
saveFlag=1
# for r in runs
#     emFibrilsToCentroids(datadir("exp_pro","cropped",r),fibrilMinSize,distanceGaussian,dilateCount,erodeCount,size1,saveFlag)
# end


