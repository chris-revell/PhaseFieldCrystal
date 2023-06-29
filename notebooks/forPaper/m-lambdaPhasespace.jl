using DrWatson; @quickactivate
using DataFrames
using CairoMakie
using Colors
using FromFile
using ImageSegmentation
using GR: delaunay
using ConcaveHull
using Images
using GeometryBasics
using LinearAlgebra
using FileIO
using ConcaveHull
using VoronoiCells
using CSV
using DataFrames
using DelimitedFiles
using LaTeXStrings
using Statistics
using Printf
using Format

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

function neighbourColours(x)
    if x == 6
        return (:white, 0.0)
    elseif x == 5
        return (:red, 1.0)
    elseif x == 7
        return (:blue, 1.0)
    else
        return (:grey, 1.0)
    end
end

function binariseSimulation!(uij)
    if uij < -0.5
        return 1.0
    else
        return 0.0
    end
end

voronoiSizeThresh = 1.3

runs = [f for f in readdir(datadir("fromCSF","lambdaTest")) if f!=".DS_Store"]

mOrder = [0.01,0.1,0.2,0.3,0.4]
λOrder = [10.0,30.0,50.0,70.0,90.0]

for mask in runs

    maskIn = load(datadir("exp_pro","masksCompressed",mask,"$mask.png"))
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

    res = collect_results(datadir("fromCSF","lambdaTest",mask))

    fig = Figure(resolution=(2000,2000),fontsize=32)

    axesDict = Dict()
    # for (i,r) in enumerate(eachrow(subset(res, :m => m -> m .== 0.1)))
    for (i,r) in enumerate(eachrow(res))
        @show i
        ax = Axis(fig[findfirst(x->x==r[:m],mOrder),findfirst(x->x==r[:λ],λOrder)],aspect=DataAspect())
        uMat = reshape(r[:u][end],(r[:nY],r[:nX]))
        
        
        
        uImg = binariseSimulation!.(uMat)
        # Segment binarised image

        seg = fast_scanning(uImg, 0.01)
        # Find centre of mass positions of all fibril segments. 
        centroidLocations = Point2{Float64}[]
        maxSize = 500
        for k in seg.segment_labels
            pixels = findall(x -> x == k, seg.image_indexmap)
            com = Tuple(sum(pixels)) ./ length(pixels)
            # Exclude the one remaining segment above size of maxSize, representing the system background
            if length(pixels) < maxSize
                push!(centroidLocations, Point2{Float64}(com[2], -com[1]))
            end
        end
        # display(map(i->maskColour(i,seg), labels_map(seg))))

        # Put centroid locations into a format for tessellation and triangulation 
        xs = [x[1] for x in centroidLocations]
        ys = [x[2] for x in centroidLocations]
        scalingFactor = maximum(size(uMat)) / (1 - 3eps(Float64))
        shiftedCentroidLocations = centroidLocations ./ scalingFactor
        shiftedCentroidLocations .+= Point2(0, 1)

        # Delaunay triangulation of centroid locations using function from GR
        n, tri = delaunay(xs, ys)
        # Count neighbours of each centroid in the triangulation 
        nNeighbours = [length(findall(x -> x == i, tri)) for i = 1:length(shiftedCentroidLocations)]

        # Concave hull to identify boundary fibrils 
        hull = concave_hull(shiftedCentroidLocations, 1)
        # Indices of fibrils within the hull 
        hullInds = sort([findall(x -> Point2(x...) == v, shiftedCentroidLocations)[1] for v in hull.vertices])

        # Voronoi tessellation of centroid positions within (0,0) (1,1) box
        tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

        tessAreas = voronoiarea(tess)
        tessAreasFiltered = [tessAreas[a] for a in 1:length(tessAreas) if a ∉ hullInds]
        # display(tessAreasFiltered)        
        meanArea = mean(tessAreasFiltered)
        # display(meanArea)

        defectCountsDict = Dict(string.(collect(3:9)) .=> zeros(Int64, 7))
        excludeCount = 0
        for j in eachindex(nNeighbours)
            if j ∉ hullInds && tessAreas[j] < voronoiSizeThresh * meanArea
                if string(nNeighbours[j]) ∈ keys(defectCountsDict)
                    defectCountsDict[string(nNeighbours[j])] += 1
                else
                    defectCountsDict[string(nNeighbours[j])] = 1
                end
            else
                excludeCount += 1
            end
        end

        runDefectProportion = 1 - defectCountsDict["6"] / (length(nNeighbours) - excludeCount) #length(hullInds))
        # push!(defectProportions, runDefectProportion)
        
        
        hidedecorations!(ax)
        hidespines!(ax)
        # Map parameters to axis in axes dictionary 
        # axesDict[(results[i, :r], results[i, :ϕ0])] = ax2

        heatmap!(ax, rotr90(uMat), colorrange=(-1.0, 1.0), colormap=:greys)
        for (j, c) in enumerate(tess.Cells)
            if j ∉ hullInds && tessAreas[j] < voronoiSizeThresh * meanArea
                vertices = [(v .- Point2(0, 1)) .* scalingFactor .+ Point2(0, size(uImg)[1]) for v in c]
                poly!(ax, vertices, color=neighbourColours(nNeighbours[j]), strokecolor=(:black, 1.0), strokewidth=1.0)
            else
                # vertices = [(v.-Point2(0,1)).*scalingFactor .+ Point2(0,size(uImg)[1]) for v in c]
                # poly!(ax2, vertices, color=:white,strokecolor=(:black,1.0),strokewidth=1.0)
            end
        end
        image!(ax, rotr90(maskImage))
        # text!(Point.([0.0], [0.0]), text=["$(@sprintf("%.3f", runDefectProportion))"], align=[(:left, :bottom)], color=:white, offset=(5, 5))
        # hullPoints = [(v.-Point2(0,1)).*scalingFactor.+Point2(0,size(uImg)[1]) for v in hull.vertices]
        # poly!(ax2,hullPoints,color=(:grey,1.0))
        # scatter!(ax2,centroidLocations.+ Point2(0,size(uImg)[1]),color=(:orange,1.0),markersize=16)
        xlims!(ax, (0, size(maskImage)[2]))
        ylims!(ax, (0, size(maskImage)[1]))
        
        Label(fig[findfirst(x->x==r[:m],mOrder),findfirst(x->x==r[:λ],λOrder),Bottom()],"m=$(r[:m]) λ=$(r[:λ]) d=$(format(runDefectProportion,precision=2))")
        
        
        
    end

    resize_to_layout!(fig)

    save(datadir("fromCSF","lambdaTest",mask,"$(mask)_lambda-m-phasespaceDefects.png"),fig)

    display(fig)

end