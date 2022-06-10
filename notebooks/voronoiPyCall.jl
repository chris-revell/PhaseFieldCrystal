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
# using GR
using GeometryBasics
using Statistics: mean
using GeometricalPredicates: Point2D
using VoronoiDelaunay


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

fileName = "$(datadir("fromCSF"))/22-06-03-21-15-58lX=244.0nX=500nY=383r=0.5tMax=250.0δt=0.01ϕ0=-0.37_finalState.png"

image = load(fileName)

grayImage = Gray.(image)

distance = 2.0

binarizedImage = binarize(grayImage,Otsu())

seg = fast_scanning(binarizedImage, 0.01)
segmentedImage = map(i->get_random_color(i), labels_map(seg))

seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<2000), (i,j)->(segment_pixel_count(seg,j)))
segmentedImage2 = map(i->get_random_color(i), labels_map(seg2))

seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)>10000), (i,j)->(segment_pixel_count(seg2,j)))
segmentedImage3 = map(i->maskColour(i,seg3), labels_map(seg3))

centroidLocations = Point2{Float64}[]
for k in seg3.segment_labels
    pixels = zeros(2)
    count = 0
    for i=1:size(image)[1]
        for j=1:size(image)[2]
            if seg3.image_indexmap[i,j] == k
                pixels .+= [j,-i]
                count += 1
            end
        end
    end
    if count<1000
        push!(centroidLocations,Point2{Float64}(pixels./count...))
    end
end

# centroidLocations .+= Point2(0.0,size(image)[1]*1.0)

# xs = [x[1] for x in centroidLocations]
# ys = [x[2] for x in centroidLocations]
# scalingFactor = maximum([xs ys])/(1-3eps(Float64))
# n, tri = delaunay(xs,ys)
# nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(centroidLocations)]

scipySpatial = pyimport("scipy.spatial")
tess = scipySpatial.Voronoi(shiftedCentroidLocations)
