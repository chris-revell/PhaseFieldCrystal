using DrWatson
@quickactivate

using Images
using ImageBinarization
using FileIO
using ImageSegmentation
using Random
using Base.Filesystem
using CairoMakie
using GR
using GeometryBasics
using LinearAlgebra
using Statistics
using ConcaveHull

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
# doubleErodeDilated = erode(erode(doubleDilate))

# Initial segmentation of eroded and dilated image
seg = fast_scanning(doubleDilate, 0.01)

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
        push!(centroidLocations,Point2(pixels./count...))
    end
end
centroidLocations .+= Point2(0.0,size(image)[1]*1.0)
xs = [x[1] for x in centroidLocations]
ys = [x[2] for x in centroidLocations]
scalingFactor = maximum([xs ys])/(1-3eps(Float64))
shiftedCentroidLocations = centroidLocations./scalingFactor

n, tri = delaunay(xs,ys)
nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

using ConcaveHull
hull = concave_hull(shiftedCentroidLocations,3)
hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])

internalCentroids = copy(shiftedCentroidLocations)
deleteat!(internalCentroids,hullInds)

pairs = Vector{Int64}[]
lengths = Float64[]
triEdges = [[1,2],[1,3],[2,3]]
for t=1:size(tri)[1]
    for e in triEdges
        if sort(tri[t,e]) ∉ pairs && tri[t,e][1] ∉ hullInds && tri[t,e][2] ∉ hullInds
            push!(pairs,sort(tri[t,e]))
            push!(lengths,norm(shiftedCentroidLocations[tri[t,e][1]].-shiftedCentroidLocations[tri[t,e][2]]))
        end
    end
end
lengths./=size(image)[2]
# deleteat!(lengths,hullInds)

using VoronoiCells
rect = Rectangle(Point2(0, 0), Point2(1, 1))
tess = voronoicells(shiftedCentroidLocations, rect)

set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(resolution=(1000,1000))

axImage = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(axImage)
image!(axImage,rotr90(image))

axSegmented = CairoMakie.Axis(fig[1,2],aspect=DataAspect())
hidedecorations!(axSegmented)
image!(axSegmented,rotr90(prunedImage2))

axTriangle = CairoMakie.Axis(fig[2,1],aspect=DataAspect())
hidedecorations!(axTriangle)
scatter!(axTriangle,shiftedCentroidLocations,color=:green,markersize=10)
xlims!(axTriangle,0,size(image)[2]/scalingFactor)
ylims!(axTriangle,0,size(image)[1]/scalingFactor)
areas = abs.(GeometryBasics.area.([shiftedCentroidLocations[tri[i,:]] for i=1:n]))
lims=(minimum(areas),maximum(areas))
for i=1:n
    poly!(axTriangle,shiftedCentroidLocations[tri[i,:]],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1.0)
end

axNeighbour = CairoMakie.Axis(fig[2,2],aspect=DataAspect())
hidedecorations!(axNeighbour)
nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]
xlims!(axNeighbour,0,size(image)[2]/scalingFactor)
ylims!(axNeighbour,0,size(image)[1]/scalingFactor)
coloursDict = Dict(2=>:black, 3=>:violet, 4=>:indigo, 5=>:blue, 6=>:white, 7=>:red, 8=>:orange)
colours = [coloursDict[x] for x in nNeighbours]
poly!(axNeighbour,hull.vertices,color=(:black,0.5))
for (i,c) in enumerate(tess.Cells)
    if i ∉ hullInds
        poly!(axNeighbour,c,color=colours[i])
    end
end

axNeighbourHistogram = CairoMakie.Axis(fig[3,1], aspect=AxisAspect(1.4))
nNeighboursFiltered = [n for (i,n) in enumerate(nNeighbours) if i ∉ hullInds]
counts = [count(x->x==y,nNeighboursFiltered) for y=minimum(nNeighbours):maximum(nNeighbours)]
barplot!(axNeighbourHistogram,collect(minimum(nNeighbours):maximum(nNeighbours)),counts)
axNeighbourHistogram.ylabel = "Count"
axNeighbourHistogram.xlabel = "Number of neighbours"

axLengths = CairoMakie.Axis(fig[3,2], aspect=AxisAspect(1.4))
hist!(axLengths,lengths,bins=50)
axLengths.xlabel = "Length"
axLengths.ylabel = "Frequency"

Label(fig[1,1,Bottom()],L"a",textsize = 48)
Label(fig[1,2,Bottom()],L"b",textsize = 48)
Label(fig[2,1,Bottom()],L"c",textsize = 48)
Label(fig[2,2,Bottom()],L"d",textsize = 48)
Label(fig[3,1,Bottom()],L"e",textsize = 48)
Label(fig[3,2,Bottom()],L"f",textsize = 48)

nVerts = length(shiftedCentroidLocations)
nEdges = length(pairs)
nFaces = length(tri)

axText = CairoMakie.Axis(fig[4,:])
text!(axText,
    "Number of fibrils = $nVerts\nNumber of edges = $nEdges\nNumber of faces = $nFaces\nMean edge length = $(mean(lengths))",
    textsize=32,
    position = Point2(0.0,0.0),
    align = (:left, :center),
    justification = :left)
hidedecorations!(axText)
hidespines!(axText)
xlims!(axText,0.0,100.0)


rowsize!(fig.layout,1,Aspect(1,1))
rowsize!(fig.layout,2,Aspect(1,1))
rowsize!(fig.layout,3,Aspect(1,1))
rowsize!(fig.layout,4,Aspect(1,0.5))
rowgap!(fig.layout,Relative(0.0))
colgap!(fig.layout,Relative(0.01))
resize_to_layout!(fig)
display(fig)
save("$(datadir("plots"))/fibrilSegmentation.png",fig)
