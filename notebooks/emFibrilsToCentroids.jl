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

function processImage(im,dilateCount)
    # Apply Gaussian filter
    distance = 1.0
    filteredImage  = imfilter(im,Kernel.gaussian(distance))
    # Binarise with Otsu algorithm
    binarizedImage = binarize(filteredImage,Otsu())
    # Erode and dilate    
    for i=1:dilateCount
        dilate!(binarizedImage)
    end
    for i=1:dilateCount
        erode!(binarizedImage)
    end
    return binarizedImage    
end

function findCentroidLocations(seg,imSize)
    centroidLocations = Point2[]
    for k in seg.segment_labels
        pixels = zeros(2)
        count = 0
        for i=1:imSize[1]
            for j=1:imSize[2]
                if seg.image_indexmap[i,j] == k
                    pixels .+= [j,-i]
                    count += 1
                end
            end
        end
        if count<1000
            push!(centroidLocations,Point2(pixels./count...))
        end
    end
    centroidLocations .+= Point2(0.0,imSize[1]*1.0)
    return centroidLocations
end

function findIntraCellSpace(dilatedImage,size1)
    # Initial segmentation of eroded and dilated image
    seg1 = fast_scanning(dilatedImage, 0.01)
    # Create inter-cellular space mask
    # Prune segments below minSizeBackground, adding pruned segments to the largest neighbouring segment     
    seg2 = prune_segments(seg1, i->(segment_pixel_count(seg1,i)<size1), (i,j)->(-segment_pixel_count(seg1,j)))
    
    intraCellSpace = (sort(collect(seg2.segment_pixel_count), by=x->x[2]))[end]
    extraCellSpace = (sort(collect(seg2.segment_pixel_count), by=x->x[2]))[end-1]
    seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<intraCellSpace.second), (i,j)->extraCellSpace.first)

    return seg1,seg2,seg3,intraCellSpace
end 

# Import image file and convert to grayscale
fileName = "/Users/christopher/Postdoc/Code/PhaseFieldCrystal/data/exp_raw/cropped/cropped_mp13ko-3wiew_4800X_hui_0002_2.png"

imageIn = load(fileName)
grayImage = Gray.(imageIn)

dilatedImage = processImage(grayImage,3)

size1 = 5000
seg1,seg2,seg3,intraCellSpace = findIntraCellSpace(dilatedImage,size1)#,size2)

newDilatedImage = copy(dilatedImage)
for i=1:size(imageIn)[1]
    for j=1:size(imageIn)[2]
        if seg3.image_indexmap[i,j] == intraCellSpace.first
            newDilatedImage[i,j] = Gray(1)
        end
    end
end

fibrilMaxSize = 1000
fibrilMinSize = 100
seg4 = fast_scanning(newDilatedImage, 0.1)
seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)>fibrilMaxSize), (i,j)->(-segment_pixel_count(seg4,j)))
seg6 = prune_segments(seg5, i->(segment_pixel_count(seg5,i)<fibrilMinSize), (i,j)->(-segment_pixel_count(seg5,j)))

centroidLocations = findCentroidLocations(seg6,size(imageIn))

fig = Figure()
ax = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
image!(ax,rotr90(imageIn))
scatter!(ax,centroidLocations,color=(:orange,0.5),markersize=10)
# save(datadir("exp_pro","segmented_$(splitpath(fileName)[end])"),fig)
display(fig)

# vals = [seg4.segment_pixel_count[i] for i in keys(seg4.segment_pixel_count)]
# segmentLabelsOrderedBySize = [i for i in keys(seg4.segment_pixel_count)]
# segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
# seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)<segment_pixel_count(seg4,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg4,j)))
# segmentedImage = map(i->maskColour(i,seg5), labels_map(seg5))
# sorted = sort(collect(seg5.segment_pixel_count), by=x->x[2])