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

function processImage(image)
    # Apply Gaussian filter
    distance = 1.0
    filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))
    # Binarise with Otsu algorithm
    binarizedImage = binarize(filteredImage,Otsu())
    # Erode and dilate
    dilatedImage = dilate(dilate(binarizedImage))
    # doubleErodeDilated = erode(erode(doubleDilate))
end


# Import image file and convert to grayscale
fileName = "/Users/christopher/Postdoc/Code/PhaseFieldCrystal/data/exp_raw/cropped/cropped_mp13ko-3wiew_4800X_hui_0002_2.png"

imageIn = load(fileName)
grayImage = Gray.(imageIn)

dilatedImage = processImage(grayImage)

# Initial segmentation of eroded and dilated image
seg1 = fast_scanning(dilatedImage, 0.01)
# Create inter-cellular space mask
# Prune segments below minSizeBackground, adding pruned segments to the largest neighbouring segment 
minSizeBackground = 1000
seg2 = prune_segments(seg1, i->(segment_pixel_count(seg1,i)<minSizeBackground), (i,j)->(-segment_pixel_count(seg1,j)))
# Prune larger segments, adding pruned segments to smallest neighbour
seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)>minSizeBackground), (i,j)->(segment_pixel_count(seg2,j)))
# vals = [seg4.segment_pixel_count[i] for i in keys(seg4.segment_pixel_count)]
# segmentLabelsOrderedBySize = [i for i in keys(seg4.segment_pixel_count)]
# segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
# seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)<segment_pixel_count(seg4,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg4,j)))
# segmentedImage = map(i->maskColour(i,seg5), labels_map(seg5))
# sorted = sort(collect(seg5.segment_pixel_count), by=x->x[2])

intraCellSpace = first.(sort(collect(seg3.segment_pixel_count), by=x->x[2]))[1]
newDilatedImage = copy(dilatedImage)
for i=1:size(imageIn)[1]
    for j=1:size(imageIn)[2]
        if seg3.image_indexmap[i,j] == intraCellSpace
            newDilatedImage[i,j] = Gray(1)
        end
    end
end

seg = fast_scanning(newDilatedImage, 0.1)
seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)>1000), (i,j)->(-segment_pixel_count(seg,j)))
seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<100), (i,j)->(-segment_pixel_count(seg2,j)))

# prunedImage2 = map(i->maskColour(i,seg3), labels_map(seg3))

centroidLocations = Point2[]
for k in seg3.segment_labels
    pixels = zeros(2)
    count = 0
    for i=1:size(imageIn)[1]
        for j=1:size(imageIn)[2]
            if seg3.image_indexmap[i,j] == k
                pixels .+= [j,-i]
                count += 1
            end
        end
    end
    if count<1000
        push!(centroidLocations,Point2(pixels./count...))
    end
end
centroidLocations .+= Point2(0.0,size(imageIn)[1]*1.0)

fig = Figure()
ax = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
image!(ax,rotr90(imageIn))
scatter!(ax,centroidLocations,color=(:orange,0.5),markersize=10)
save(datadir("exp_pro","segmented_$(splitpath(fileName)[end])"),fig)
display(fig)

