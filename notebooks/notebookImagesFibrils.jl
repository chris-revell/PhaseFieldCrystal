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

seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<30000), (i,j)->(-segment_pixel_count(seg,j)))
segmentedImage .= map(i->maskColour(i,seg2), labels_map(seg2))
