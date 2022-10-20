using DrWatson
@quickactivate
using UnPack
using CairoMakie; set_theme!(figure_padding=0,fontsize=32)
using ColorSchemes
using DifferentialEquations
using JLD2
using LaTeXStrings
using DataFrames
using Images
using ImageBinarization
using ImageSegmentation
using ConcaveHull
using VoronoiCells
using GR: delaunay

# Location of data within data/fromCSF/ 
folderPath = "PhaseSpace3"

# Collate results as a dataframe 
results = collect_results(datadir("fromCSF",folderPath); subfolders = true)

# Set up figure canvas and dictionary to map parameters to axes
fig = Figure(resolution=(6000,6000))
axesDict = Dict()

coloursDict = Dict(1=>:black,2=>:black, 3=>:violet, 4=>:indigo, 5=>:blue, 6=>:white, 7=>:red, 8=>:orange, 9=>:orange, 10=>:orange, 11=>:orange, 12=>:orange, 13=>:orange, 14=>:orange)

# Set max size for fibril segments 
maxSize = 500

# Loop to process data from each run 
for i=1:nrow(subset(results, :m => m -> m.== 0.1))
    display(i)
    # Convert simulation result to a 2D matrix
    uMat = reshape(results[i,:u],(results[i,:nY],results[i,:nX]))#transpose(reshape(results[i,:u],(results[i,:nY],results[i,:nX])))
    # Convert matrix to a grayscale image
    uImg = Gray.(uMat)
    # Binarise grayscale image
    uBinarized = binarize(uImg,Otsu())
    # Segment binarised image
    seg = fast_scanning(uBinarized, 0.01)    
    # Prune segments larger than maxsize, adding their pixels to the largest neighbouring segment.
    # Should leave segments for each fibril and one segment for the background.
    seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)>500), (i,j)->(-segment_pixel_count(seg,j)))
    
    # Find centre of mass positions of all fibril segments. 
    centroidLocations = Point2{Float64}[]
    for k in seg2.segment_labels
        pixels = zeros(2)
        count = 0
        for i=1:size(uImg)[1]
            for j=1:size(uImg)[2]
                if seg2.image_indexmap[i,j] == k
                    pixels .+= [j,-i]
                    count += 1
                end
            end
        end
        # Exclude the one remaining segment above size of maxSize, representing the system background
        if count<maxSize
            push!(centroidLocations,Point2{Float64}(pixels./count...))
        end
    end
    # Rescale centroid locations into a (0,0) (1,1) box for voronoi tessellation
    centroidLocations .+= Point2(0.0,size(uImg)[1]*1.0)
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum([xs ys])/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor

    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs,ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

    # Concave hull to identify boundary fibrils 
    hull = concave_hull(shiftedCentroidLocations,3)
    # Indices of fibrils within the hull 
    hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])

    # Voronoi tessellation of centroid positions within (0,0) (1,1) box
    tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

    # Map neighbour counts to colours 
    colours = [coloursDict[x] for x in nNeighbours]

    # New axis within figure
    ax = CairoMakie.Axis(fig[i%6+1,i÷6+1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    # Map parameters to axis in axes dictionary 
    axesDict[(results[i,:r],results[i,:ϕ0])] = ax
    
    # Plot centroid locations and triangulation 
    # scatter!(ax,shiftedCentroidLocations,color=:green,markersize=10)
    # for i=1:n
    #     poly!(ax,shiftedCentroidLocations[tri[i,:]],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1.0)
    # end
    
    # Plot voronoi cells around each fibril coloured by neighbour count
    poly!(ax,hull.vertices,color=(:black,0.5))
    for (i,c) in enumerate(tess.Cells)
        if i ∉ hullInds
            poly!(ax,c,color=colours[i])
        end
    end
    xlims!(ax,0,size(uImg)[2]/scalingFactor)
    ylims!(ax,0,size(uImg)[1]/scalingFactor)
end 

# Rearrange axes to put relevant r and ϕ0 values in order
sortedrs  = unique!(sort(first.(keys(axes))))
sortedϕ0s = unique!(sort(last.(keys(axes))))
for (i,r) in enumerate(sortedrs)
    for (j,ϕ0) in enumerate(sortedϕ0s)
        fig[length(sortedrs)+1-i,length(sortedϕ0s)+1-j][1,1] = axesDict[(r,ϕ0)]
        # Colorbar(fig[length(sortedrs)-i,length(sortedϕ0s)-j][1,2], limits=(-1,1),colormap=:bwr)#, size = 25)
        Label(fig[length(sortedrs)+1-i,length(sortedϕ0s)+1-j,Bottom()],L"r=%$r, \phi_0=%$ϕ0",textsize=64)
    end
end

# Resize columns
for i=1:length(sortedϕ0s)
    colsize!(fig.layout,i,Aspect(1.0,1.0))
end
resize_to_layout!(fig)

display(fig)

save(datadir("fromCSF",folderPath,"neighbourPhaseSpace.png"),fig)


# internalCentroids = deleteat!(shiftedCentroidLocations,hullInds) # copy(shiftedCentroidLocations)
    # pairs = Vector{Int64}[]
    # lengths = Float64[]
    # triEdges = [[1,2],[1,3],[2,3]]
    # for t=1:size(tri)[1]
    #     for e in triEdges
    #         if sort(tri[t,e]) ∉ pairs
    #             if !(tri[t,e][1] ∈ hullInds && tri[t,e][2] ∈ hullInds)
    #                 push!(pairs,sort(tri[t,e]))
    #                 push!(lengths,norm(shiftedCentroidLocations[tri[t,e][1]].-shiftedCentroidLocations[tri[t,e][2]]))
    #             end
    #         end
    #     end
    # end
    # lengths./=size(image)[2]
    # deleteat!(lengths,hullInds)