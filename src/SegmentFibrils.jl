module SegmentFibrils

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
using VoronoiCells

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

function segmentFibrils(fileName,distance,nDilate,nErode,thresh1,thresh2,thresh3,displayFlag,saveFlag)

    # Import image file and convert to grayscale
    image = load(fileName)
    grayImage = Gray.(image)

    # Apply Gaussian filter

    filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))

    # Binarise with Otsu algorithm
    binarizedImage = binarize(filteredImage,Otsu())

    # Erode and dilate
    for i=1:nDilate
        dilate!(binarizedImage)
    end
    for i=1:nErode
        erode!(binarizedImage)
    end

    # Initial segmentation of eroded and dilated image
    seg = fast_scanning(binarizedImage, 0.1)

    # Prune segments below
    seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<thresh1), (i,j)->(segment_pixel_count(seg,j)))
    vals = [seg2.segment_pixel_count[i] for i in keys(seg2.segment_pixel_count)]
    segmentLabelsOrderedBySize = [i for i in keys(seg2.segment_pixel_count)]
    segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
    seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<segment_pixel_count(seg2,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg2,j)))
    segmentedImage = map(i->maskColour(i,seg3), labels_map(seg3))
    sorted = sort(collect(seg3.segment_pixel_count), by=x->x[2])

    newGrayImage = binarizedImage
    for i=1:size(image)[1]
        for j=1:size(image)[2]
            if seg3.image_indexmap[i,j] == sorted[1][1]
                newGrayImage[i,j] = Gray(1)
            end
        end
    end

    seg = fast_scanning(newGrayImage, 0.1)
    seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)>thresh2), (i,j)->(-segment_pixel_count(seg,j)))
    seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<thresh3), (i,j)->(-segment_pixel_count(seg2,j)))

    prunedImage = map(i->maskColour(i,seg3), labels_map(seg3))

    centroidLocations = Point2[]
    for k in seg3.segment_labels
        pixels = zeros(2)
        count = 0
        for i=1:size(prunedImage)[1]
            for j=1:size(prunedImage)[2]
                if seg3.image_indexmap[i,j] == k
                    pixels .+= [j,-i]
                    count += 1
                end
            end
        end
        # Don't include the centroid of the largest (background) segment
        if count<1000
            push!(centroidLocations,Point2(pixels./count...))
        end
    end
    # Rescale centroid locations to be between (0,0) and (1,1)
    centroidLocations .+= Point2(0.0,size(image)[1]*1.0)
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum([xs ys])/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor

    # Perform delaunay triangulation to return number of triangles and the set of indices for each triangle
    n, tri = delaunay(xs,ys)
    # Find number of neighbours for each fibril
    nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

    # Find concave hull using 3 neighbours to identify boundary fibrils
    hull = concave_hull(shiftedCentroidLocations,3)
    # Hack to find indices of boundary fibrils from locations
    hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])

    # Exclude boundary fibril locations
    internalCentroids = copy(shiftedCentroidLocations)
    deleteat!(internalCentroids,hullInds)

    # Find neighbouring fibril pairs and corresponding lengths
    pairs = Vector{Int64}[]
    lengths = Float64[]
    triEdges = [[1,2],[1,3],[2,3]]  # All possible neighbour pairs in a triangle
    boundaryEdges = Int64[]
    pairsCount = 0
    for t=1:size(tri)[1]
        for e in triEdges
            # Skip if this pair of edges already exists in the pairs list
            if sort(tri[t,e]) ∉ pairs
                pairsCount += 1
                push!(pairs,sort(tri[t,e]))
                push!(lengths,norm(shiftedCentroidLocations[tri[t,e][1]].-shiftedCentroidLocations[tri[t,e][2]]))
                # Make a note if both of the fibrils for this pair are in the concave hull
                if tri[t,e][1] ∈ hullInds && tri[t,e][2] ∈ hullInds
                    push!(boundaryEdges,pairsCount)
                end
            end
        end
    end
    lengths./=size(image)[2]
    # Create another vector of lengths, removing those pairs for which both fibrils are at the periphery
    lengthsFiltered = copy(lengths)
    deleteat!(lengthsFiltered,boundaryEdges)

    # Perform voronoi tesselation of fibril locations
    rect = Rectangle(Point2(0, 0), Point2(1, 1))
    tess = voronoicells(shiftedCentroidLocations, rect)

    # Create figure canvas
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig = Figure(resolution=(1000,1000))

    # Display original image in first axis
    axImage = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(axImage)
    image!(axImage,rotr90(image))

    # Display segmented image in second axis
    axSegmented = CairoMakie.Axis(fig[1,2],aspect=DataAspect())
    hidedecorations!(axSegmented)
    image!(axSegmented,rotr90(prunedImage))

    # Display delaunay triangulation in third axis
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

    # Display voronoi tesselation with voronoi cells coloured by fibril neighbour count in 4th axis. Exclude fibrils that are part of concave hull.
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

    # Display histogram of fibril neighbour count in 5th axis, excluding fibrils that are part of the concave hull.
    axNeighbourHistogram = CairoMakie.Axis(fig[3,1], aspect=AxisAspect(1.4))
    nNeighboursFiltered = [n for (i,n) in enumerate(nNeighbours) if i ∉ hullInds]
    counts = [count(x->x==y,nNeighboursFiltered) for y=minimum(nNeighbours):maximum(nNeighbours)]
    barplot!(axNeighbourHistogram,collect(minimum(nNeighbours):maximum(nNeighbours)),counts)
    axNeighbourHistogram.ylabel = "Count"
    axNeighbourHistogram.xlabel = "Number of neighbours"

    # Display histogram of fibril neighbour pair lengths in 6th axis, excluding fibrils that are part of the concave hull.
    axLengths = CairoMakie.Axis(fig[3,2], aspect=AxisAspect(1.4))
    hist!(axLengths,lengthsFiltered,bins=50)
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
        "Number of fibrils = $nVerts\nNumber of edges = $nEdges\nNumber of faces = $nFaces\nMean edge length = $(mean(lengthsFiltered))",
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

    displayFlag==1 ? display(fig) : nothing

    saveFlag ==1 ? save("$(datadir("plots"))/fibrilSegmentation_$fileName",fig) : nothing
end

export segmentFibrils

end
