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
using CairoMakie
using GR
using GeometryBasics


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

filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))

binarizedImage = binarize(filteredImage,Otsu())

doubleDilate = dilate(dilate(dilate(binarizedImage)))

# doubleErodeDilated = erode(erode(doubleDilate))

seg = fast_scanning(doubleErodeDilated, 0.01)

seg4 = prune_segments(seg, i->(segment_pixel_count(seg,i)<1000), (i,j)->(segment_pixel_count(seg,j)))
vals = [seg4.segment_pixel_count[i] for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize = [i for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)<segment_pixel_count(seg4,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg4,j)))
segmentedImage = map(i->get_random_color(i), labels_map(seg5))
sorted = sort(collect(seg5.segment_pixel_count), by=x->x[2])

newGrayImage = doubleDilate
for i=1:size(image)[1]
    for j=1:size(image)[2]
        if seg5.image_indexmap[i,j] == sorted[1][1]
            newGrayImage[i,j] = Gray(1)
        end
    end
end

seg = fast_scanning(newGrayImage, 0.1)
seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)>1000), (i,j)->(-segment_pixel_count(seg,j)))
seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<100), (i,j)->(-segment_pixel_count(seg2,j)))

# segmentedImage = map(i->get_random_color(i), labels_map(seg))
# prunedImage1 = map(i->maskColour(i,seg2), labels_map(seg2))
# prunedImage2 = map(i->maskColour(i,seg3), labels_map(seg3))
# segmentLocations = Dict(seg3.segment_labels .=> fill(Tuple[],length(seg3.segment_labels)))
# for i=1:size(prunedImage2)[1]
#     for j=1:size(prunedImage2)[2]
#         segmentLocations[seg3.image_indexmap[i,j]] = push!(segmentLocations[seg3.image_indexmap[i,j]],(i,j))
#         # display(seg3.image_indexmap[i,j])
#     end
# end


centroidLocations = Point2[]
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
        push!(centroidLocations,Point2(pixels./count))
    end
end

fig = Figure(); ax = CairoMakie.Axis(fig[1,1],aspect=DataAspect())

scatter!(ax,centroidLocations)

xs = [x[1] for x in centroidLocations]
ys = [x[2] for x in centroidLocations]
n, tri = delaunay(xs,ys)

for i=1:n
    poly!(ax,centroidLocations[tri[i,:]],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1.0)
end

display(fig)

using LazySets
hull = convex_hull(Vector.(centroidLocations))
