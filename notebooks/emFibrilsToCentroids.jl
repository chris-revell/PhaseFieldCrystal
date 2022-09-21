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
@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions


function emFibrilsToCentroids(fileName,fibrilMinSize,distance,dilateCount,erodeCount)
    mkpath(datadir("exp_pro","emCentroids"))
    # Import image file and convert to grayscale
    imageIn = load(fileName)
    imSize = size(imageIn)
    grayImage = Gray.(imageIn)
    # Apply Gaussian filter
    filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))
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

    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    image!(ax,rotr90(imageIn))
    scatter!(ax,centroidLocations,color=(:orange,1.0),markersize=6)
    save(datadir("exp_pro","emCentroids","segmented_$(splitpath(fileName)[end])"),fig)
    save(datadir("exp_pro","emCentroids","segmented_$(splitpath(fileName)[end][1:end-4]).jld2"),@strdict fileName fibrilMinSize distance dilateCount erodeCount centroidLocations fig)
    display(fig)
    return fig, centroidLocations, newIndexMap
end

fileName = "/Users/christopher/Postdoc/Code/PhaseFieldCrystal/data/exp_raw/cropped/cropped_mp13ko-3wiew_4800X_hui_0002_2.png"
distance = 1.0
fibrilMinSize = 250
dilateCount = 2
erodeCount = 2
emFibrilsToCentroids(fileName,fibrilMinSize,distance,dilateCount,erodeCount)


# function processImage(im,dilateCount,distance)
#     # Apply Gaussian filter
#     filteredImage  = imfilter(im,Kernel.gaussian(distance))
#     # Binarise with Otsu algorithm
#     binarizedImage = binarize(filteredImage,Otsu())
#     # Erode and dilate    
#     for i=1:dilateCount
#         dilate!(binarizedImage)
#     end
#     for i=1:dilateCount
#         erode!(binarizedImage)
#     end
#     return binarizedImage    
# end

# function findCentroidLocations(indexMap,imSize,extraCellSpace)
#     # Find segment labels within index map
#     labels = unique(indexMap)
#     # Remove 0 segment, used to represent intra-cell space
#     filter!(val->val≠0,labels)
#     # Remove extraCellSpace segment, used to represent extra-cellular space
#     filter!(val->val≠extraCellSpace,labels)
#     centroidLocations = Point2[]
#     for k in labels
#         pixels = zeros(2)
#         count = 0
#         for i=1:imSize[1]
#             for j=1:imSize[2]
#                 if newIndexMap[i,j] == k
#                     pixels .+= [j,-i]
#                     count += 1
#                 end
#             end
#         end
#         count>250 ? push!(centroidLocations,Point2(pixels./count)) : nothing
#     end
#     centroidLocations .+= Point2(0.0,imSize[1]*1.0)
#     return centroidLocations
# end

# function findIntraCellSpace(dilatedImage,size1)
#     # Initial segmentation of eroded and dilated image
#     seg1 = fast_scanning(dilatedImage, 0.01)
#     # Create inter-cellular space mask
#     # Prune segments below minSizeBackground, adding pruned segments to the largest neighbouring segment     
#     seg2 = prune_segments(seg1, i->(segment_pixel_count(seg1,i)<size1), (i,j)->(-segment_pixel_count(seg1,j)))
#     seg3 = prune_segments(seg2, first.(sort(collect(seg2.segment_pixel_count), by=x->x[2])[1:end-2]), (i,j)->-segment_pixel_count(seg2,j))

#     return seg1,seg2,seg3
# end 




# vals = [seg4.segment_pixel_count[i] for i in keys(seg4.segment_pixel_count)]
# segmentLabelsOrderedBySize = [i for i in keys(seg4.segment_pixel_count)]
# segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
# seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)<segment_pixel_count(seg4,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg4,j)))
# segmentedImage = map(i->maskColour(i,seg5), labels_map(seg5))
# sorted = sort(collect(seg5.segment_pixel_count), by=x->x[2])