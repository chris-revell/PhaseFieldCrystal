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

function binariseSimulation!(uij)
    if uij > 0.5
        return 1.0
    else
        return 0.0
    end
end 

results = collect_results!(datadir("fromCSF","allMasks"); subfolders = true)

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

fig = Figure(resolution=(6000,6000),fontsize=64)

axes = Dict()
sizes = Dict()

for (i,r) in enumerate(runs)
    
    subsetResults = subset(results, :path => m -> occursin.(r[1:end-4],m))

    maskIn = load(datadir("exp_pro","masksCompressed",r[1:end-4],"$(r[1:end-4]).png"))
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
    uMat = reshape(subsetResults[1,:u][end],(subsetResults[1,:nY],subsetResults[1,:nX]))
    # Binarise grayscale image
    uImg = binariseSimulation!.(uMat)
    # Convert matrix to a grayscale image
    uGray = Gray.(uImg)
    # Segment binarised image
    seg = fast_scanning(uGray, 0.01)    
    # Find centre of mass positions of all fibril segments. 
    centroidLocations = Point2{Float64}[]
    maxSize = 500
    for k in seg.segment_labels
        pixels = findall(x->x==k,seg.image_indexmap)
        com = Tuple(sum(pixels))./length(pixels)        
        # Exclude the one remaining segment above size of maxSize, representing the system background
        if length(pixels)<maxSize
            push!(centroidLocations,Point2{Float64}(com[2],-com[1]))
        end
    end
    # display(map(i->get_random_color(i), seg.image_indexmap))

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(abs.([xs ys]))/(1-3eps(Float64))
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
    tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

    ax2 = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect())

    for (i,c) in enumerate(tess.Cells)
        if i ∉ hullInds
            vertices = [(v.-Point2(0,1)).*scalingFactor .+ Point2(0,size(uImg)[1]) for v in c]
            poly!(ax2, vertices, color=neighbourColours(nNeighbours[i]),strokecolor=(:black,1.0),strokewidth=1.0)
        else
            # poly!(ax2, c.*scalingFactor, color=:white,strokecolor=(:black,1.0),strokewidth=1.0)
        end
    end
    image!(ax2,rotr90(maskImage))
    poly!(ax2,hull.vertices,color=(:grey,1.0))
    scatter!(ax2,centroidLocations.+ Point2(0,size(uImg)[1]),color=(:orange,1.0),markersize=4)
    hidedecorations!(ax2)
    hidespines!(ax2)
    
    Label(fig[(i-1)%6+1,(i-1)÷6+1, BottomLeft()], "$i", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0), color=:black)
    
    axes[r] = ax2
    sizes[r] = size(maskImage)
end 

xMax = maximum(first.(values(sizes)))
yMax = maximum(last.(values(sizes)))

for r in runs 
    xlims!(axes[r],(0,xMax))
    ylims!(axes[r],(0,yMax))
end

colgap!(fig.layout, 1, -700)
colgap!(fig.layout, 2, -600)
colgap!(fig.layout, 3, -600)
colgap!(fig.layout, 4, -500)
colgap!(fig.layout, 5, -100)

rowgap!(fig.layout, 1, -100)
rowgap!(fig.layout, 2, -100)
rowgap!(fig.layout, 3, -200)

resize_to_layout!(fig)
display(fig)

save(datadir("fromCSF","allMasks","allMasksSimulationNeighbourGrid.png"),fig)
