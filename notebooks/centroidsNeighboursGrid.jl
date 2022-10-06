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

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

fig = Figure(resolution=(6000,6000))

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


mkpath(datadir("exp_pro","emCentroidNeighbours"))

for (i,r) in enumerate(runs)
    maskData = load(datadir("exp_pro","masks",r[1:end-4],"$(r[1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData

    centroidData = load(datadir("exp_pro","emCentroidsInteractive","$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum([xs ys])/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor

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

    maskImage = fill(RGBA(1,1,1,1),size(newIndexMap))
    for i=1:size(newIndexMap)[1]
        for j=1:size(newIndexMap)[2]
            if newIndexMap[i,j] < 0.5
                maskImage[i,j] = RGBA(0.0,0.0,0.0,1.0)
            else
                maskImage[i,j] = RGBA(0.0,0.0,0.0,0.0)
            end
        end
    end    

    ax2 = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect())

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
end 


resize_to_layout!(fig)
display(fig)

save(datadir("exp_pro","emCentroidNeighbours","grid.png"),fig)


