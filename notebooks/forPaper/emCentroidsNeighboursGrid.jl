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
using DelimitedFiles
using Printf
using Statistics
using StatsBase
using DelimitedFiles

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

mkpath(datadir("exp_pro", "emCentroidNeighbours"))

runs = [r for r in Vector(readdlm(datadir("exp_pro", "filesToUse.txt"))[:, 1]) if !(occursin("mp13ko", r) || occursin("18tailT_4800X_HUI_0007_0", r) || occursin("18tailT_4800X_HUI_0008_0", r))]

croppedLX = DataFrame(CSV.File(datadir("exp_pro", "lengthMeasurements", "croppedLX.csv")))

voronoiSizeThresh = 1.5

fig = Figure(resolution=(6000, 6000), fontsize=64)

runsToUse = [1, 2, 12, 13, 17, 18]

defectProportions = Float64[]

defectCountsDataFrame = DataFrame()

for (k, r) in enumerate(runs)#[runsToUse])

    maskData = load(datadir("exp_pro", "masks", r[1:end-4], "$(r[1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData

    centroidData = load(datadir("exp_pro", "emCentroidsInteractive", "$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData

    imageIn = load(datadir("exp_pro", "cropped", r))
    imSize = size(imageIn)
    grayImage = Gray.(imageIn)

    croppedLXRow = filter(:file => f -> f == r, croppedLX)
    lengthPerPixel = croppedLXRow[1, :lX] / size(imageIn)[2]

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(size(newIndexMap)) / (1 - 3eps(Float64))
    shiftedCentroidLocations = centroidLocations ./ scalingFactor

    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs, ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x -> x == i, tri)) for i = 1:length(shiftedCentroidLocations)]

    # Concave hull to identify boundary fibrils 
    hull = concave_hull(shiftedCentroidLocations)
    # Indices of fibrils within the hull 
    hullInds = sort([findall(x -> Point2(x...) == v, shiftedCentroidLocations)[1] for v in hull.vertices])

    # Voronoi tessellation of centroid positions within (0,0) (1,1) box
    tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

    tessAreas = voronoiarea(tess)
    tessAreasFiltered = [tessAreas[a] for a in 1:length(tessAreas) if a ∉ hullInds]
    # display(tessAreasFiltered)        
    meanArea = mean(tessAreasFiltered)
    # display(meanArea)

    defectCountsDict = Dict(string.(collect(4:8)) .=> zeros(Int64, 5))
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
    runDefectProportion = 1 - defectCountsDict["6"] / (length(nNeighbours) - excludeCount)
    push!(defectProportions, runDefectProportion)

    dfLocal = DataFrame(defectCountsDict)
    dfLocal[!,:file] = [r]
    dfLocal[!,:k] = [k]
    dfLocal[!,:defectProportion] = [runDefectProportion]
    defectCountsDataFrame = vcat(defectCountsDataFrame,dfLocal,cols = :union)

    ax = CairoMakie.Axis(fig[(k-1)÷6+1, mod(k-1,6)+1], aspect=DataAspect(), backgroundcolor=:white)

    image!(ax, rotr90(grayImage))

    for (j, c) in enumerate(tess.Cells)
        if j ∉ hullInds && tessAreas[j] < voronoiSizeThresh * meanArea
            vertices = [v .* scalingFactor for v in c]
            poly!(ax, vertices, color=neighbourColours(nNeighbours[j]), strokecolor=(:black, 1.0), strokewidth=1.0)
        end
    end

    if k == 26 || k == 27
        text!(Point.([0.0], [0.0]), text=["$(@sprintf("%.3f", runDefectProportion))"], align=[(:left, :bottom)], color=:white, offset=(5, 5), fontsize=64)
    else
        text!(Point.([0.0], [size(grayImage)[1]]), text=["$(@sprintf("%.3f", runDefectProportion))"], align=[(:left, :top)], color=:white, offset=(5, 5), fontsize=64)
    end

    hidedecorations!(ax)
    hidespines!(ax)

    Label(fig[(k-1)÷6+1, mod(k-1,6)+1, Bottom()], "$k", padding=(0, 10, 10, 0), color=:black, fontsize=128)

end

CSV.write(datadir("exp_pro", "emCentroidMeasurements", "defectCounts.csv"), defectCountsDataFrame)

resize_to_layout!(fig)
display(fig)

save(datadir("exp_pro", "emCentroidNeighbours", "emNeighboursGridWithDefectProportionUpdatedSubset.png"), fig)

# points = Point2[]
# hNorm = normalize(fit(Histogram, defectProportions, 0:0.1:1))
# edgeVec = collect(hNorm.edges[1])
# for i = 1:length(hNorm.weights)
#     push!(points, Point2(mean(edgeVec[i:i+1]), hNorm.weights[i]))
# end
fig = CairoMakie.Figure(resolution=(500,500), fontsize=32)
ax = CairoMakie.Axis(fig[1,1])
hist!(ax, simDefectPropotions; bins=collect(0:0.1:1), normalization=:none, color=(:blue,0.75), label="Simulation")
hist!(ax,defectProportions; bins=collect(0:0.1:1), normalization=:none, color = (:green,0.75), label="EM")
ax.xlabel = "Defect proportion"
ax.ylabel = "Frequency"
axislegend(ax, position=:rt, labelsize=18)
display(fig)
save(datadir("exp_pro", "emCentroidNeighbours", "emAndSimDefectProportionsHistogram.png"), fig)

writedlm(datadir("exp_pro", "emCentroidNeighbours","defectProportions.txt"),defectProportions)