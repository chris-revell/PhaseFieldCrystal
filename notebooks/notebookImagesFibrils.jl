using DrWatson
@quickactivate

using Images
using ImageBinarization
using FileIO
using ImageSegmentation
using ImageTransformations
using ImageView
using ImageSegmentation
using Random
using Base.Filesystem
using ImageSmooth



function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

function maskColour(i,seg)
    if seg.segment_pixel_count[i]==maximum(values(seg.segment_pixel_count))
        return RGB{}(0,0,0)
    else
        return RGB{}(1,1,1)
    end
end

fileName = "data/exp_raw/Cropped1_mmp13ko-3wiew_4800X_hui_0002.png"

image = load(fileName)

grayImage = Gray.(image)

distance = 2.0

filteredImage  = imfilter(Gray.(image),Kernel.gaussian(distance))

binarizedImage = binarize(grayImage,Otsu())

seg = fast_scanning(binarizedImage, 0.1)

biggestSegment = sort(collect(seg.segment_pixel_count), by=x->x[2])[end][1]

seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)>1200), (i,j)->(biggestSegment))
seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<500), (i,j)->(biggestSegment))
segmentedImage = map(i->maskColour(i,seg), labels_map(seg))
prunedImage1 = map(i->maskColour(i,seg2), labels_map(seg2))
prunedImage2 = map(i->maskColour(i,seg3), labels_map(seg3))
