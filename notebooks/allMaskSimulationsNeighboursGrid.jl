using DrWatson; @quickactivate
using Images
using CairoMakie
using Colors
using GeometryBasics
using LinearAlgebra
using FromFile
using FileIO
using ConcaveHull
using VoronoiCells
using GR: delaunay
using CSV
using DataFrames
using ImageBinarization
using ImageSegmentation

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

results = collect_results!(datadir("fromCSF","allMasks"); subfolders = true)

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

fig = Figure(resolution=(6000,6000),fontsize=64)

axes = Dict()
sizes = Dict()

fig = Figure(resolution=(6000,6000),backgroundcolor=:black,fontsize=64)

function neighbourColours(x)
    if x==6
        return :white
    elseif x==5 
        return :red 
    elseif x==7
        return :blue
    else 
        return :black 
    end
end

for (i,r) in enumerate(runs)
    
    subsetResults = subset(results, :path => m -> occursin.(r[1:end-4],m))

    maskIn = load(datadir("exp_pro","masksCompressed",r[1:end-4],"$(r[1:end-4]).png"))
    # @unpack newIndexMap, lX, h = maskData
    maskImage = fill(RGBA(1,1,1,1),size(maskIn))
    for i=1:size(maskIn)[1]
        for j=1:size(maskIn)[2]
            if maskIn[i,j] < 0.5
                maskImage[i,j] = RGBA(0.0,0.0,0.0,1.0)
            else
                maskImage[i,j] = RGBA(0.0,0.0,0.0,0.0)
            end
        end
    end    

    # Convert simulation result to a 2D matrix
    uMat = reshape(subsetResults[1,:u][end],(subsetResults[1,:nY],subsetResults[1,:nX]))#transpose(reshape(results[i,:u],(results[i,:nY],results[i,:nX])))
    # Convert matrix to a grayscale image
    uImg = Gray.(uMat)
    # Binarise grayscale image
    uBinarized = binarize(uImg,Otsu())
    # Segment binarised image
    seg = fast_scanning(uBinarized, 0.01)    
    # Prune segments larger than maxsize, adding their pixels to the largest neighbouring segment.
    # Should leave segments for each fibril and one segment for the background.
    seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)>500), (i,j)->(-segment_pixel_count(seg,j)))
    seg2 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<50), (i,j)->(-segment_pixel_count(seg2,j)))
    display(map(i->get_random_color(i), labels_map(seg2)))
    # Find centre of mass positions of all fibril segments. 
    centroidLocations = Point2{Float64}[]
    maxSize = 500
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



    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum([xs ys])/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor
    shiftedCentroidLocations .+= Point2(0,1)

    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs,ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

    # Concave hull to identify boundary fibrils 
    hull = concave_hull(shiftedCentroidLocations,1)
    # Indices of fibrils within the hull 
    hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])

    # Voronoi tessellation of centroid positions within (0,0) (1,1) box
    # display(shiftedCentroidLocations)
    tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

    ax2 = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect(), backgroundcolor=:black)

    for (i,c) in enumerate(tess.Cells)
        if i ∉ hullInds
            poly!(ax2, c.*scalingFactor, color=neighbourColours(nNeighbours[i]),strokecolor=(:black,1.0),strokewidth=1.0)
        else
            # poly!(ax2, c.*scalingFactor, color=:white,strokecolor=(:black,1.0),strokewidth=1.0)
        end
    end
    image!(ax2,rotr90(maskImage))
    # poly!(ax2,hull.vertices,color=(:grey,1.0))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=4)
    hidedecorations!(ax2)
    hidespines!(ax2)
    
    Label(fig[(i-1)%6+1,(i-1)÷6+1, Bottom()], "$i", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0), color=:white)
    
    axes[r] = ax2
    sizes[r] = size(maskImage)
end 


lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))

lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
lengthPerPixelDict = Dict()
lengthDict = Dict()
for r in runs 
    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
    lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
    lengthDict[r] = lengthPerPixelDict[r].*sizes[r]
end

xMax = maximum(first.(values(lengthDict)))
yMax = maximum(last.(values(lengthDict)))

for r in runs 
    xlims!(axes[r],(0,xMax)./lengthPerPixelDict[r])
    ylims!(axes[r],(0,yMax)./lengthPerPixelDict[r])
end
# for r in runs 
#     xlims!(axes[r],(-xMax,xMax)./(2*lengthPerPixelDict[r]).+first(sizes[r])/2)
#     ylims!(axes[r],(-yMax,yMax)./(2*lengthPerPixelDict[r]).+last(sizes[r])/2)
# end


resize_to_layout!(fig)
display(fig)

save(datadir("exp_pro","emCentroidNeighbours","grid.png"),fig)



