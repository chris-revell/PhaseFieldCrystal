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

function neighbourColours(x)
    if x == 6
        return (:white, 0.0)
    elseif x == 5
        return (:red, 0.75)
    elseif x == 7
        return (:blue, 0.75)
    else
        return (:grey, 0.75)
    end
end

function binariseSimulation!(uij)
    if uij < -0.5
        return 1.0
    else
        return 0.0
    end
end

# fig = Figure(fontsize=32,resolution=(800,500))
# ax = Axis(fig[1,1])
# ylims!(ax,(0,1))
r = "17tailT_4800X_HUI_0002_0"
voronoiSizeThresh = 1.3

runID = 5

maskIn = load(datadir("exp_pro", "masksCompressed", r, "$r.png"))
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

results = collect_results(datadir("sims", "timeResolution", r))

nX = results[runID, :nX]
nY = results[runID, :nY]

timepointsShort = collect(1:15:80)#[1, 6, 11, 16, 21, 26]
timepointsLong = collect(101:180:1001)#[101, 281, 461, 551, 701, 801]

fig1 = Figure(figure_padding=0, resolution=(6000, 2000), fontsize=64)

for (i, t) in enumerate(timepointsShort)

    uMat = reshape(results[runID, :u][t], (nY, nX))

    # fig1 = Figure(figure_padding=0,resolution=(1000,2000),fontsize=64)
    ax1 = CairoMakie.Axis(fig1[1, i], aspect=DataAspect())
    uInternal = Observable(zeros(nX, nY))
    heatmap!(ax1, transpose(uMat), colorrange=(-1.0, 1.0), colormap=:bwr)
    hidedecorations!(ax1)
    hidespines!(ax1)
    ax1.yreversed = true
    image!(ax1, transpose(maskImage))
    ax1.title = "t=$(@sprintf("%.2f", results[runID,:t][t]))"

    # fig3 = Figure(figure_padding=0,resolution=(1000,1000),fontsize=64)
    ax2 = CairoMakie.Axis(fig1[2, i], aspect=DataAspect())
    heatmap!(ax2, rotr90(uMat), colorrange=(-1.0, 1.0), colormap=:greys)
    hidedecorations!(ax2)
    hidespines!(ax2)
    # ax2.title = "t=$(@sprintf("%.2f", results[runID,:t][i]))"
    image!(ax2, rotr90(maskImage))
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
    scalingFactor = maximum(size(uInternal[])) / (1 - 3eps(Float64))
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
        if "6" in keys(defectCountsDict)
            runDefectProportion = 1 - defectCountsDict["6"] / (length(nNeighbours) - excludeCount)
        else
            runDefectProportion = 1
        end
        for (j, c) in enumerate(tess.Cells)
            if j ∉ hullInds && tessAreas[j] < voronoiSizeThresh * meanArea
                vertices = [(v .- Point2(0, 1)) .* scalingFactor .+ Point2(0, size(uImg)[1]) for v in c]
                poly!(ax2, vertices, color=neighbourColours(nNeighbours[j]), strokecolor=(:black, 1.0), strokewidth=1.0)
            end
        end
        text!(Point.([0.0], [0.0]), text=["$(@sprintf("%.3f", runDefectProportion))"], align=[(:left, :bottom)], color=:white, offset=(5, 5))
    end
    xlims!(ax2, (0, size(maskImage)[2]))
    ylims!(ax2, (0, size(maskImage)[1]))
    # resize_to_layout!(fig1)
    # save(datadir("sims","timeResolution",r,"timeResolutionBoth$(i).png"),fig1)
end
resize_to_layout!(fig1)
save(datadir("sims", "timeResolution", r, "timeResolutionShort$runID.png"), fig1)

fig1 = Figure(figure_padding=0, resolution=(6000, 2000), fontsize=64)

for (i, t) in enumerate(timepointsLong)

    uMat = reshape(results[runID, :u][t], (nY, nX))

    # fig1 = Figure(figure_padding=0,resolution=(1000,2000),fontsize=64)
    ax1 = CairoMakie.Axis(fig1[1, i], aspect=DataAspect())
    uInternal = Observable(zeros(nX, nY))
    heatmap!(ax1, transpose(uMat), colorrange=(-1.0, 1.0), colormap=:bwr)
    hidedecorations!(ax1)
    hidespines!(ax1)
    ax1.yreversed = true
    image!(ax1, transpose(maskImage))
    ax1.title = "t=$(@sprintf("%.2f", results[runID,:t][t]))"

    # fig3 = Figure(figure_padding=0,resolution=(1000,1000),fontsize=64)
    ax2 = CairoMakie.Axis(fig1[2, i], aspect=DataAspect())
    heatmap!(ax2, rotr90(uMat), colorrange=(-1.0, 1.0), colormap=:greys)
    hidedecorations!(ax2)
    hidespines!(ax2)
    # ax2.title = "t=$(@sprintf("%.2f", results[runID,:t][i]))"
    image!(ax2, rotr90(maskImage))
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
    scalingFactor = maximum(size(uInternal[])) / (1 - 3eps(Float64))
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
            if j ∉ hullInds && tessAreas[j] < 1.3 * meanArea
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
        for (j, c) in enumerate(tess.Cells)
            if j ∉ hullInds && tessAreas[j] < 1.3 * meanArea
                vertices = [(v .- Point2(0, 1)) .* scalingFactor .+ Point2(0, size(uImg)[1]) for v in c]
                poly!(ax2, vertices, color=neighbourColours(nNeighbours[j]), strokecolor=(:black, 1.0), strokewidth=1.0)
            end
        end
        text!(Point.([0.0], [0.0]), text=["$(@sprintf("%.3f", runDefectProportion))"], align=[(:left, :bottom)], color=:white, offset=(5, 5))
    end
    xlims!(ax2, (0, size(maskImage)[2]))
    ylims!(ax2, (0, size(maskImage)[1]))
end
resize_to_layout!(fig1)
save(datadir("sims", "timeResolution", r, "timeResolutionLong$runID.png"), fig1)
