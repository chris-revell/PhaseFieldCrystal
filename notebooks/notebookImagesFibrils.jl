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
using GeometryBasics
using CairoMakie

function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

function maskColour(i,seg)
    if seg.segment_pixel_count[i]==maximum(values(seg.segment_pixel_count))
        return Gray{}(0)
    else
        return Gray{}(1)
    end
end

fileName = "data/exp_raw/Cropped1_mmp13ko-3wiew_4800X_hui_0002.png"

image = load(fileName)

grayImage = Gray.(image)

distance = 1.0

filteredImage  = imfilter(Gray.(image),Kernel.gaussian(distance))

binarizedImage = binarize(filteredImage,Otsu())

# closedImage = closing(binarizedImage)

doubleDilate = dilate(dilate(dilate(binarizedImage)))

doubleErodeDilated = erode(erode(doubleDilate))

# mgImage = morpholaplace(binarizedImage)

# seg = felzenszwalb(doubleDilate, 300,100)
seg = fast_scanning(doubleDilate, 0.01)

# biggestSegment = sort(collect(seg.segment_pixel_count), by=x->x[2])[end][1]
# diff_fun(i,neighbour) = -segment_pixel_count(seg,neighbour)

seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)>1000), (i,j)->(-segment_pixel_count(seg,j)))
seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<100), (i,j)->(-segment_pixel_count(seg2,j)))
segmentedImage = map(i->Gray(values(seg.segment_pixel_count[i])/maximum(values(seg.segment_pixel_count))), labels_map(seg))
# prunedImage1 = map(i->maskColour(i,seg2), labels_map(seg2))
prunedImage2 = map(i->maskColour(i,seg3), labels_map(seg3))


seg4 = prune_segments(seg, i->(segment_pixel_count(seg,i)<1000), (i,j)->(-segment_pixel_count(seg,j)))
seg5 = prune_segments(seg, i->(segment_pixel_count(seg,i)<39000), (i,j)->(-segment_pixel_count(seg,j)))
segmentedImage5 = map(i->get_random_color(i), labels_map(seg5))


# segmentLocations = Dict(seg3.segment_labels .=> fill(Tuple[],length(seg3.segment_labels)))
# for i=1:size(prunedImage2)[1]
#     for j=1:size(prunedImage2)[2]
#         segmentLocations[seg3.image_indexmap[i,j]] = push!(segmentLocations[seg3.image_indexmap[i,j]],(i,j))
#         # display(seg3.image_indexmap[i,j])
#     end
# end

centroidLocations = Point2f[]
for k in seg3.segment_labels
    pixels = zeros(2)
    count = 0
    for i=1:size(prunedImage2)[1]
        for j=1:size(prunedImage2)[2]
            if seg3.image_indexmap[i,j] == k
                pixels .+= [j,-i]
                count += 1
            end
        end
    end
    if count<1000
        push!(centroidLocations,Point2f(pixels./count))
    end
end

fig = CairoMakie.Figure(); ax = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
scatter!(ax,centroidLocations)

display(fig)
