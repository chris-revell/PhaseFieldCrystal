using DrWatson;
@quickactivate;
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
using DelimitedFiles
using Statistics
using Printf

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

function filterFunction(r, ϕ0)
    r == 0.8 && ϕ0 == 0.4
end

# runs = [r for r in Vector(readdlm(datadir("exp_pro", "filesToUse.txt"))[:, 1]) if !(occursin("mp13ko", r) || occursin("18tailT_4800X_HUI_0007_0", r) || occursin("18tailT_4800X_HUI_0008_0", r))]

runs = [r for r in Vector(readdlm(datadir("exp_pro", "filesToUse.txt"))[:, 1]) if !(occursin("mp13ko", r) || occursin("18tailT_4800X_HUI_0007_0", r) || occursin("18tailT_4800X_HUI_0008_0", r))]

voronoiSizeThresh = 1.3

fig = Figure(resolution=(6000, 6000), fontsize=64)

for (i, r) in enumerate(runs)

    results = collect_results(datadir("fromCSF", "allMasksPhasespaceSeparateLengths", r[1:end-4]); subfolders=true)

    subsetResults = filter([:r, :ϕ0] => filterFunction, results) #subset(results, :path => m -> occursin.(r[1:end-4],m))

    maskIn = load(datadir("exp_pro", "masksCompressed", r[1:end-4], "$(r[1:end-4]).png"))
    maskImage = fill(RGBA(1, 1, 1, 1), size(maskIn))
    for i = 1:size(maskIn)[1]
        for j = 1:size(maskIn)[2]
            if maskIn[i, j] < 0.5
                maskImage[i, j] = RGBA(0.0, 0.0, 0.0, 1.0)
            else
                maskImage[i, j] = RGBA(0.0, 0.0, 0.0, 0.0)
            end
        end
    end

    # Convert simulation result to a 2D matrix
    uMat = reshape(subsetResults[1, :u][end], (subsetResults[1, :nY], subsetResults[1, :nX]))
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
        pixels = findall(x -> x == k, seg.image_indexmap)
        com = Tuple(sum(pixels)) ./ length(pixels)
        # Exclude the one remaining segment above size of maxSize, representing the system background
        if length(pixels) < maxSize
            push!(centroidLocations, Point2{Float64}(com[2], -com[1]))
        end
    end
    # display(map(i->get_random_color(i), seg.image_indexmap))

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(abs.([xs ys])) / (1 - 3eps(Float64))
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

    defectCountsDict = Dict(string.(collect(3:9)).=>zeros(Int64,7))
    excludeCount = 0
    for j in eachindex(nNeighbours)
        if j ∉ hullInds && tessAreas[j] < voronoiSizeThresh*meanArea
            if string(nNeighbours[j]) ∈ keys(defectCountsDict)
                defectCountsDict[string(nNeighbours[j])] += 1
            else
                defectCountsDict[string(nNeighbours[j])] = 1
            end
        else
            excludeCount += 1
        end
    end
    runDefectProportion = 1 - defectCountsDict["6"]/(length(nNeighbours)-excludeCount)

    ax = CairoMakie.Axis(fig[(i-1)÷6+1, (i-1)%6+1], aspect=DataAspect())

    heatmap!(ax, rotr90(uMat), colorrange=(-1.0, 1.0), colormap=:grays)
    for (j, c) in enumerate(tess.Cells)
        if j ∉ hullInds && tessAreas[j] < voronoiSizeThresh*meanArea
            vertices = [(v .- Point2(0, 1)) .* scalingFactor .+ Point2(0, size(uImg)[1]) for v in c]
            poly!(ax, vertices, color=neighbourColours(nNeighbours[j]), strokecolor=(:black, 1.0), strokewidth=1.0)
        end
    end
    image!(ax, rotr90(maskImage))
    if i == 26 || i == 27
        text!(Point.([0.0], [0.0]), text=["$(@sprintf("%.3f", runDefectProportion))"], align=[(:left, :bottom)], color=:white, offset=(5, 5), fontsize=64)
    else
        text!(Point.([0.0], [size(maskImage)[1]]), text=["$(@sprintf("%.3f", runDefectProportion))"], align=[(:left, :top)], color=:white, offset=(5, 5), fontsize=64)
    end
    hidedecorations!(ax)
    hidespines!(ax)
    xlims!(ax, (0, size(maskImage)[2]))
    ylims!(ax, (0, size(maskImage)[1]))

    Label(fig[(i-1)÷6+1, (i-1)%6+1, Bottom()], "$i", padding=(0, 10, 10, 0), color=:black, fontsize=128)

end

resize_to_layout!(fig)
display(fig)

save(datadir("fromCSF", "allMasksPhasespaceSeparateLengths", "allMasksSimulationNeighbourGrid.png"), fig)
