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
EMfileName = "data/exp_raw/Cropped1_mmp13ko-3wiew_4800X_hui_0002.png"
EMimage = load(EMfileName)
EMgrayImage = Gray.(EMimage)

# Apply Gaussian filter
distance = 1.0
EMfilteredImage  = imfilter(EMgrayImage,Kernel.gaussian(distance))

# Binarise with Otsu algorithm
EMbinarizedImage = binarize(EMfilteredImage,Otsu())

# Erode and dilate
EMdoubleDilate = dilate(dilate(dilate(EMbinarizedImage)))
# doubleErodeDilated = erode(erode(doubleDilate))

# Initial segmentation of eroded and dilated image
EMseg = fast_scanning(EMdoubleDilate, 0.01)

# Prune segments below
EMseg4 = prune_segments(EMseg, i->(segment_pixel_count(EMseg,i)<1000), (i,j)->(segment_pixel_count(EMseg,j)))
EMvals = [EMseg4.segment_pixel_count[i] for i in keys(EMseg4.segment_pixel_count)]
EMsegmentLabelsOrderedBySize = [i for i in keys(EMseg4.segment_pixel_count)]
EMsegmentLabelsOrderedBySize .= EMsegmentLabelsOrderedBySize[sortperm(EMvals)]
EMseg5 = prune_segments(EMseg4, i->(segment_pixel_count(EMseg4,i)<segment_pixel_count(EMseg4,EMsegmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(EMseg4,j)))
EMsegmentedImage = map(i->maskColour(i,EMseg5), labels_map(EMseg5))
EMsorted = sort(collect(EMseg5.segment_pixel_count), by=x->x[2])

EMnewGrayImage = EMdoubleDilate
for i=1:size(EMimage)[1]
    for j=1:size(EMimage)[2]
        if EMseg5.image_indexmap[i,j] == EMsorted[1][1]
            EMnewGrayImage[i,j] = Gray(1)
        end
    end
end

EMseg = fast_scanning(EMnewGrayImage, 0.1)
EMseg2 = prune_segments(EMseg, i->(segment_pixel_count(EMseg,i)>1000), (i,j)->(-segment_pixel_count(EMseg,j)))
EMseg3 = prune_segments(EMseg2, i->(segment_pixel_count(EMseg2,i)<100), (i,j)->(-segment_pixel_count(EMseg2,j)))

EMprunedImage2 = map(i->maskColour(i,EMseg3), labels_map(EMseg3))

EMcentroidLocations = Point2[]
for k in EMseg3.segment_labels
    pixels = zeros(2)
    count = 0
    for i=1:size(EMprunedImage2)[1]
        for j=1:size(EMprunedImage2)[2]
            if EMseg3.image_indexmap[i,j] == k
                pixels .+= [j,-i]
                count += 1
            end
        end
    end
    if count<1000
        push!(EMcentroidLocations,Point2(pixels./count...))
    end
end
EMcentroidLocations .+= Point2(0.0,size(EMimage)[1]*1.0)
EMxs = [x[1] for x in EMcentroidLocations]
EMys = [x[2] for x in EMcentroidLocations]
EMscalingFactor = maximum([EMxs EMys])/(1-3eps(Float64))
EMshiftedCentroidLocations = EMcentroidLocations./EMscalingFactor

EMn, EMtri = delaunay(EMxs,EMys)
EMnNeighbours = [length(findall(x->x==i,EMtri)) for i=1:length(EMshiftedCentroidLocations)]

using ConcaveHull
EMhull = concave_hull(EMshiftedCentroidLocations,3)
EMhullInds = sort([findall(x->Point2(x...)==v,EMshiftedCentroidLocations)[1] for v in EMhull.vertices])

EMinternalCentroids = copy(EMshiftedCentroidLocations)
deleteat!(EMinternalCentroids,EMhullInds)

EMpairs = Vector{Int64}[]
EMlengths = Float64[]
triEdges = [[1,2],[1,3],[2,3]]
for t=1:size(EMtri)[1]
    for e in triEdges
        push!(EMpairs,sort(EMtri[t,e]))
        if !(EMtri[t,e[1]] ∈ EMhullInds || EMtri[t,e[2]] ∈ EMhullInds)
            push!(EMlengths,norm(EMshiftedCentroidLocations[EMtri[t,e][1]].-EMshiftedCentroidLocations[EMtri[t,e][2]]))
        end
    end
end
EMlengths./=size(EMimage)[2]
# deleteat!(lengths,hullInds)

using VoronoiCells
EMrect = Rectangle(Point2(0, 0), Point2(1, 1))
EMtess = voronoicells(EMshiftedCentroidLocations, EMrect)





SimfileName = "$(datadir("fromCSF"))/22-06-09-19-14-13lX=200.0nX=500nY=383r=0.5tMax=500.0δt=0.01ϕ0=-0.37_finalState.png"
Simimage = load(SimfileName)
SimgrayImage = Gray.(Simimage)
distance = 2.0
SimbinarizedImage = binarize(SimgrayImage,Otsu())
Simseg = fast_scanning(SimbinarizedImage, 0.01)
SimsegmentedImage = map(i->get_random_color(i), labels_map(Simseg))
Simseg2 = prune_segments(Simseg, i->(segment_pixel_count(Simseg,i)<2000), (i,j)->(segment_pixel_count(Simseg,j)))
SimsegmentedImage2 = map(i->get_random_color(i), labels_map(Simseg2))
Simseg3 = prune_segments(Simseg2, i->(segment_pixel_count(Simseg2,i)>10000), (i,j)->(segment_pixel_count(Simseg2,j)))
SimsegmentedImage3 = map(i->maskColour(i,Simseg3), labels_map(Simseg3))
SimcentroidLocations = Point2{Float64}[]
for k in Simseg3.segment_labels
    pixels = zeros(2)
    count = 0
    for i=1:size(Simimage)[1]
        for j=1:size(Simimage)[2]
            if Simseg3.image_indexmap[i,j] == k
                pixels .+= [j,-i]
                count += 1
            end
        end
    end
    if count<1000
        push!(SimcentroidLocations,Point2{Float64}(pixels./count...))
    end
end
SimcentroidLocations .+= Point2(0.0,size(Simimage)[1]*1.0)
Simxs = [x[1] for x in SimcentroidLocations]
Simys = [x[2] for x in SimcentroidLocations]
SimscalingFactor = maximum([Simxs Simys])/(1-3eps(Float64))
SimshiftedCentroidLocations = SimcentroidLocations./SimscalingFactor

Simn, Simtri = delaunay(Simxs,Simys)
SimnNeighbours = [length(findall(x->x==i,Simtri)) for i=1:length(SimshiftedCentroidLocations)]

using ConcaveHull
Simhull = concave_hull(SimshiftedCentroidLocations,3)
SimhullInds = sort([findall(x->Point2(x...)==v,SimshiftedCentroidLocations)[1] for v in Simhull.vertices])

SiminternalCentroids = copy(SimshiftedCentroidLocations)
deleteat!(SiminternalCentroids,SimhullInds)

Simpairs = Vector{Int64}[]
Simlengths = Float64[]
triEdges = [[1,2],[1,3],[2,3]]
for t=1:size(Simtri)[1]
    for e in triEdges
        if sort(Simtri[t,e]) ∉ Simpairs
            push!(Simpairs,sort(Simtri[t,e]))
            if !(Simtri[t,e[1]] ∈ SimhullInds || Simtri[t,e[2]] ∈ SimhullInds)
                push!(Simlengths,norm(SimshiftedCentroidLocations[Simtri[t,e][1]].-SimshiftedCentroidLocations[Simtri[t,e][2]]))
            end
        end
    end
end
Simlengths./=size(Simimage)[2]
# deleteat!(lengths,hullInds)
Simrect = Rectangle(Point2(0, 0), Point2(1, 1))
Simtess = voronoicells(SimshiftedCentroidLocations, Simrect)


set_theme!(figure_padding=4, backgroundcolor=(:white,1.0), font="Helvetica")
# fig = Figure()
# axImage = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
# hidedecorations!(axImage)
# image!(axImage,rotr90(Simimage))
# resize_to_layout!(fig)
#

# axSegmented = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
# hidedecorations!(axSegmented)
# image!(axSegmented,rotr90(SimsegmentedImage3))


# axTriangle = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
# hidedecorations!(axTriangle)
# scatter!(axTriangle,EMshiftedCentroidLocations,color=:green,markersize=10)
# xlims!(axTriangle,0,size(EMimage)[2]/EMscalingFactor)
# ylims!(axTriangle,0,size(EMimage)[1]/EMscalingFactor)
# # areas = abs.(area.([ shiftedCentroidLocations[tri[i,:]] for i=1:n]))
# # lims=(minimum(areas),maximum(areas))
# for i=1:EMn
#     poly!(axTriangle,EMshiftedCentroidLocations[EMtri[i,:]],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1.0)
# end


# axNeighbour = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
# hidedecorations!(axNeighbour)
# nNeighbours = [length(findall(x->x==i,Simtri)) for i=1:length(SimshiftedCentroidLocations)]
# xlims!(axNeighbour,0,size(Simimage)[2]/SimscalingFactor)
# ylims!(axNeighbour,0,size(Simimage)[1]/SimscalingFactor)
# coloursDict = Dict(2=>:black, 3=>:violet, 4=>:indigo, 5=>:blue, 6=>:white, 7=>:red, 8=>:orange)
# colours = [coloursDict[x] for x in nNeighbours]
# poly!(axNeighbour,Simhull.vertices,color=(:black,0.5))
# for (i,c) in enumerate(Simtess.Cells)
#     if i ∉ SimhullInds
#         poly!(axNeighbour,c,color=colours[i])
#     end
# end
fig = Figure()
SimnNeighbours = [length(findall(x->x==i,Simtri)) for i=1:length(SimshiftedCentroidLocations)]
axNeighbourHistogram = CairoMakie.Axis(fig[1,1], aspect=AxisAspect(1.4))
SimnNeighboursFiltered = [n for (i,n) in enumerate(SimnNeighbours) if i ∉ SimhullInds]
Simcounts = [count(x->x==y,SimnNeighboursFiltered) for y=minimum(SimnNeighbours):maximum(SimnNeighbours)]

EMnNeighbours = [length(findall(x->x==i,EMtri)) for i=1:length(EMshiftedCentroidLocations)]
EMnNeighboursFiltered = [n for (i,n) in enumerate(EMnNeighbours) if i ∉ EMhullInds]
EMcounts = [count(x->x==y,EMnNeighboursFiltered) for y=minimum(EMnNeighbours):maximum(EMnNeighbours)]
barplot!(axNeighbourHistogram,collect(minimum(EMnNeighbours):maximum(EMnNeighbours)).+0.2,EMcounts,label="Electron micrograph",width=0.5)
barplot!(axNeighbourHistogram,collect(minimum(SimnNeighbours):maximum(SimnNeighbours)).-0.2,Simcounts,label="Simulation",width=0.5)
axNeighbourHistogram.ylabel = "Count"
axNeighbourHistogram.xlabel = "Number of neighbours"
xlims!(axNeighbourHistogram,3,9)
axislegend(axNeighbourHistogram)
# axLengths = CairoMakie.Axis(fig[3,2], aspect=AxisAspect(1.4))
# hist!(axLengths,lengths,bins=50)
# axLengths.xlabel = "Length"
# axLengths.ylabel = "Frequency"


save("$(projectdir())/EMvsSimnNeighbours.png",fig)

#
# fig = Figure()
# axLengths = CairoMakie.Axis(fig[1,1], aspect=AxisAspect(1.4))
# hist!(axLengths,EMlengths,bins=50,label="Electron micrograph")
# hist!(axLengths,Simlengths,bins=50,label="Simulation")
# axLengths.xlabel = "Length"
# axLengths.ylabel = "Frequency"
# xlims!(axLengths,2.0E-5,4.0E-5)
# axislegend(axLengths)
#
#
# display(fig)
# save("$(projectdir())/EMvsSimLengths.png",fig)
