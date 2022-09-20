using DrWatson; @quickactivate
using FromFile
using Images
using ImageBinarization
using FileIO
using ImageSegmentation
using ImageSegmentation
using Random
using Base.Filesystem
using CairoMakie
using GR
using GeometryBasics
using ImageView
@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

fileName = "data/exp_pro/cropped/cropped_mp13ko-3wiew_4800X_hui_0002_0.png"

imageIn = load(fileName)

grayImage = Gray.(imageIn)

distance = 1.0

filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))

binarizedImage = binarize(filteredImage,Intermodes())

dilate!(binarizedImage)
dilate!(binarizedImage)
dilate!(binarizedImage)

erode!(binarizedImage)
erode!(binarizedImage)
erode!(binarizedImage)

seg = fast_scanning(binarizedImage, 0.01)

segmentedImage = map(i->get_random_color(i), labels_map(seg))

seg4 = prune_segments(seg, i->(segment_pixel_count(seg,i)<1000), (i,j)->(segment_pixel_count(seg,j)))
vals = [seg4.segment_pixel_count[i] for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize = [i for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)<segment_pixel_count(seg4,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg4,j)))
prunedImage = map(i->maskColour(i,seg5), labels_map(seg5))
ImageView.imshow(prunedImage)

# save(datadir("exp_pro","masks",splitpath(fileName)[end],prunedImage) #Need to resize too

new_width = 500
percentage_scale = new_width/size(imageIn,2)
new_size = trunc.(Int64, size(prunedImage) .* percentage_scale)
img_rescaled = imresize(prunedImage, new_size)
save(datadir("exp_pro","masks",splitpath(fileName)[end]),img_rescaled)

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
