using DrWatson;
@quickactivate;
using DataFrames
using CairoMakie
using Colors
using FromFile
using GeometryBasics
using NumericalIntegration
using Printf
using GR: delaunay
using VoronoiCells
using ImageSegmentation
using ConcaveHull
using Statistics

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

function binariseSimulation!(uij)
    if uij < -0.5
        return 1.0
    else
        return 0.0
    end
end

r = "17tailT_4800X_HUI_0002_0"
runID = 6

voronoiSizeThresh = 1.3

results = collect_results(datadir("sims", "timeResolution", r))
nY = results[runID, :nY]
nX = results[runID, :nX]

points = Point2[]

for i = 1:5:1001#:5:101

    uMat = reshape(results[runID, :u][i], (nY, nX))
    # Binarise grayscale image
    uImg = binariseSimulation!.(uMat)
    seg = fast_scanning(uImg, 0.01)
    centroidLocations = Point2{Float64}[]
    maxSize = 500
    for k in seg.segment_labels
        pixels = findall(x -> x == k, seg.image_indexmap)
        com = Tuple(sum(pixels)) ./ length(pixels)
        if length(pixels) < maxSize
            push!(centroidLocations, Point2{Float64}(com[2], -com[1]))
        end
    end
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(size(uMat)) / (1 - 3eps(Float64))
    shiftedCentroidLocations = centroidLocations ./ scalingFactor
    shiftedCentroidLocations .+= Point2(0, 1)
    if length(shiftedCentroidLocations) > 2
        n, tri = delaunay(xs, ys)
        nNeighbours = [length(findall(x -> x == i, tri)) for i = 1:length(shiftedCentroidLocations)]
        hull = concave_hull(shiftedCentroidLocations, 1)
        hullInds = sort([findall(x -> Point2(x...) == v, shiftedCentroidLocations)[1] for v in hull.vertices])
        tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))
        tessAreas = voronoiarea(tess)
        tessAreasFiltered = [tessAreas[a] for a in 1:length(tessAreas) if a ∉ hullInds]
        meanArea = mean(tessAreasFiltered)
        defectCountsDict = Dict()
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
        runDefectProportion = 1 - defectCountsDict["6"] / (length(nNeighbours) - excludeCount)
        push!(points, Point2(results[runID, :t][i], runDefectProportion))
    end
end

fig1 = Figure(figure_padding=30, resolution=(500, 500), fontsize=32)
ax1 = CairoMakie.Axis(fig1[1, 1],xscale=Makie.log10)
lines!(ax1, points)
ax1.xlabel = "Time"
ax1.ylabel = "Defect proportion"
xlims!(ax1,(1,1000))
colsize!(fig1.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig1)
display(fig1)
save(datadir("sims", "timeResolution", r, "defectsAgainstTimeLog$(runID).png"), fig1)