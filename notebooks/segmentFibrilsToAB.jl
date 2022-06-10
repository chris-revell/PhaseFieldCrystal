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
using SparseArrays

# Function to set random colour for each segment
function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

# Function to set white for the largest segment and black for others
function maskColour(i,seg)
    if seg.segment_pixel_count[i]==maximum(values(seg.segment_pixel_count))
        return Gray{}(1)
    else
        return Gray{}(0)
    end
end

# Import image file and convert to grayscale
fileName = "data/exp_raw/Cropped1_mmp13ko-3wiew_4800X_hui_0002.png"
image = load(fileName)
grayImage = Gray.(image)

# Apply Gaussian filter
distance = 1.0
filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))

# Binarise with Otsu algorithm
binarizedImage = binarize(filteredImage,Otsu())

# Erode and dilate
doubleDilate = dilate(dilate(dilate(binarizedImage)))
doubleErodeDilated = erode(erode(doubleDilate))

# Initial segmentation of eroded and dilated image
seg = fast_scanning(doubleErodeDilated, 0.01)

# Prune segments below
seg4 = prune_segments(seg, i->(segment_pixel_count(seg,i)<1000), (i,j)->(segment_pixel_count(seg,j)))
vals = [seg4.segment_pixel_count[i] for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize = [i for i in keys(seg4.segment_pixel_count)]
segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
seg5 = prune_segments(seg4, i->(segment_pixel_count(seg4,i)<segment_pixel_count(seg4,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg4,j)))
segmentedImage = map(i->maskColour(i,seg5), labels_map(seg5))
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

prunedImage2 = map(i->maskColour(i,seg3), labels_map(seg3))

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

xs = [x[1] for x in centroidLocations]
ys = [x[2] for x in centroidLocations]
n, tri = delaunay(xs,ys)
nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(centroidLocations)]

pairs = Vector{Int64}[]
# lengths = Float64[]
triEdges = [[1,2],[1,3],[2,3]]
for t=1:size(tri)[1]
    for e in triEdges
        if sort(tri[t,e]) ∉ pairs
            push!(pairs,sort(tri[t,e]))
            # push!(lengths,norm(centroidLocations[tri[t,e][1]].-centroidLocations[tri[t,e][2]]))
        end
    end
end
# lengths./=size(image)[2]


nVerts = length(centroidLocations)
nEdges = length(pairs)
nFaces = size(tri)[1]

A = spzeros(Int64,nEdges,nVerts)
B = spzeros(Int64,nFaces,nEdges)

for (i,p) in enumerate(pairs)
    display(i)
    display(p)
    A[i,p[1]] = -1
    A[i,p[2]] = 1
end
