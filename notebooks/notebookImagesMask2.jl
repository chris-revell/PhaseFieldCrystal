using DrWatson
@quickactivate

using Images
using ImageBinarization
using FileIO
using ImageSegmentation
using ImageTransformations
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
        return Gray{}(1)
    else
        return Gray{}(0)
    end
end

fileName = "data/exp_raw/Cropped1_mmp13ko-3wiew_4800X_hui_0002.png"

image = load(fileName)

grayImage = Gray.(image)

distance = 1.0

filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))

binarizedImage = binarize(filteredImage,Otsu())

doubleDilate = dilate(dilate(dilate(binarizedImage)))

doubleErodeDilated = erode(erode(doubleDilate))

seg = fast_scanning(doubleErodeDilated, 0.01)

segmentedImage = map(i->get_random_color(i), labels_map(seg))

seg4 = prune_segments(seg, i->(segment_pixel_count(seg,i)<1000), (i,j)->(segment_pixel_count(seg,j)))
vals = [seg4.segment_pixel_count[i] for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize = [i for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)<segment_pixel_count(seg4,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg4,j)))
prunedImage = map(i->maskColour(i,seg5), labels_map(seg5))

save(datadir("exp_pro","mask_$(Filesystem.splitpath(fileName)[end])"),prunedImage) #Need to resize too

# set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica",fontsize=48)
# fig = Figure(resolution=(1000,1000))
# axImage = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
# hidedecorations!(axImage)
# hidespines!(axImage)
# image!(axImage,rotr90(image))
# axSegmented = CairoMakie.Axis(fig[1,2],aspect=DataAspect())
# hidedecorations!(axSegmented)
# hidespines!(axSegmented)
# image!(axSegmented,rotr90(segmentedImage))
# axPruned = CairoMakie.Axis(fig[1,3],aspect=DataAspect())
# hidedecorations!(axPruned)
# hidespines!(axPruned)
# image!(axPruned,rotr90(prunedImage))
# Label(fig[1,1,Bottom()],L"a",textsize = 48)
# Label(fig[1,2,Bottom()],L"b",textsize = 48)
# Label(fig[1,3,Bottom()],L"c",textsize = 48)
# rowsize!(fig.layout,1,Aspect(1,1))
# resize_to_layout!(fig)
# display(fig)
# save("maskSegmentation.png",fig)
